source("workflow/scripts/helpers/find_neighbors.R")
source("workflow/scripts/helpers/connect_KNN_igraph.R")

options(scipen=999)

library(FNN)
library(ComplexHeatmap)

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

input_path         <- args$input_path
input_path_rse     <- args$input_path_rse
input_path_counts  <- args$input_path_counts
train_samples_path <- args$train_samples
C                  <- as.integer(args$C)
j                  <- as.numeric(args$j)
output_path_plot   <- args$output_path_plot
output_path_map    <- args$output_path_map

TRAIN <- readLines(train_samples_path)

dir.create(dirname(output_path_map),  recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_path_plot), recursive = TRUE, showWarnings = FALSE)

parts <- unlist(strsplit(basename(input_path), split = "_"))
TxID <- parts[length(parts)]

# CHANGE: go straight to the abs matrix used for clustering and drop the
# data.frame immediately. We no longer keep dataset / dataset_clustering /
# DC as three separate n_cpg x n_cpg objects at once -- just DC.
dataset <- read.delim(input_path, check.names = FALSE, header = TRUE, sep = "\t")
DC <- abs(as.matrix(dataset))
rm(dataset); invisible(gc(verbose = FALSE))

CpGIDs <- colnames(DC)
parts_split <- strsplit(CpGIDs, "_", fixed = TRUE)
chrs   <- vapply(parts_split, `[`, character(1), 1)
starts <- vapply(parts_split, `[`, character(1), 2)
coords <- data.frame(chr = chrs, start = starts, end = starts, CpGID = CpGIDs)

# ---- load expression for this TxID (TRAIN counts) ----
counts <- read.table(input_path_counts, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
expr   <- as.numeric(counts[TxID, ])
names(expr) <- colnames(counts)
rm(counts); invisible(gc(verbose = FALSE))

# ---- load methylation matrix from rse for this TxID ----
# CHANGE: compute common_samples against the rse's colnames (lazy, doesn't
# touch data) and subset the SummarizedExperiment BEFORE materializing to a
# dense matrix, so as.matrix() only ever reads the TRAIN columns off HDF5
# instead of loading every sample and immediately discarding most of them.
rse <- HDF5Array::loadHDF5SummarizedExperiment(input_path_rse)

# LEAKAGE-SAFE: the j-term is scored on TRAIN samples only.
# mean_abs_cor_expr selects the clustering parameters, so it must not see
# test methylation or test expression.
common_samples <- Reduce(intersect, list(colnames(rse), names(expr), TRAIN))
cat(sprintf("[%s] j-term scored on %d TRAIN samples\n", TxID, length(common_samples)))
if (j != 0 && length(common_samples) < 3)
  stop(sprintf("[%s] only %d train samples for the j-term.", TxID, length(common_samples)))

meth <- as.matrix(SummarizedExperiment::assay(rse[, common_samples], "M"))
rm(rse); invisible(gc(verbose = FALSE))
expr <- expr[common_samples]
invisible(gc(verbose = FALSE))

# =============================================================================
# Clustering grid search
#
# Same computation as the original nested loops, reorganised for speed and
# memory:
#   - find_neighbors() depends only on (coords, i)  -> computed ONCE per CpG
#   - quantile()/diff() depend only on (CpG, percentile) -> all flex_points are
#     read off the SAME sorted diff vector
#   - each grid point is SCORED IMMEDIATELY and its igraph object discarded,
#     so only the cluster lists survive
# =============================================================================
n_cpg  <- nrow(coords)

percentiles_list <- seq(from = 10, to = 1000, length.out = 100)
flex_points_list <- seq(from = 1, to = 3, by = 1)

cat(sprintf("[%s] %d CpGs, %d grid points\n",
            TxID, n_cpg, length(percentiles_list) * length(flex_points_list)))

neighbors.spatial.list <- lapply(seq_len(n_cpg), function(i) find_neighbors(coords, i, 1, TRUE))
row_medians            <- matrixStats::rowMedians(DC)

# ---- presort each row of DC ONCE ----------------------------------------
# Previously: quantile(DC[i, ], probs) re-sorts row i from scratch, once per
# percentile_denom -> up to 100 redundant sorts of the same row. And
# which(DC[i, ] >= threshold) does a fresh O(n_cpg) scan for every one of the
# 300 grid points. Both only ever need the row's sort order, which is fixed.
# Sorting all rows once here and reusing that turns:
#   - threshold computation from O(100 * n_cpg^2 log n_cpg) into O(n_cpg^2 log n_cpg)
#   - neighbor lookup from O(300 * n_cpg^2) into O(n_cpg^2 log n_cpg) + O(300 * n_cpg log n_cpg)
# Trade-off: this holds two extra n_cpg x n_cpg-sized objects (order indices
# + sorted values) alongside DC for the duration of the grid search --
# roughly another ~1.5x DC's memory footprint. If that's tight for your
# largest transcripts, drop C for this rule rather than skip this section.
row_order_desc <- vector("list", n_cpg)   # indices, DC value descending
row_sorted_asc <- vector("list", n_cpg)   # DC values, ascending (for quantile)
for (i in seq_len(n_cpg)) {
    ord <- order(DC[i, ], decreasing = TRUE)
    row_order_desc[[i]] <- ord
    row_sorted_asc[[i]] <- rev(DC[i, ord])
}

# stats::quantile(x, probs) with default type-7 interpolation, computed
# directly off an already-sorted vector -- no re-sort.
quantile_from_sorted <- function(sorted_asc, probs) {
    n  <- length(sorted_asc)
    h  <- (n - 1) * probs + 1
    lo <- pmin(pmax(floor(h), 1), n)
    hi <- pmin(pmax(ceiling(h), 1), n)
    sorted_asc[lo] + (h - lo) * (sorted_asc[hi] - sorted_asc[lo])
}

score_grid_point <- function(clusters, graph, membership, key) {
    n_clusters    <- length(unique(membership))
    modularity    <- igraph::modularity(graph, membership)
    if (!is.finite(modularity)) modularity <- 0      # no edges -> undefined
    modularity    <- round(modularity, 3)
    fragmentation <- round(n_clusters / igraph::vcount(graph), 3)

    cluster_id_of <- rep(seq_along(clusters), times = lengths(clusters))
    cpg_of_cluster <- coords$CpGID[unlist(clusters)]

    keep <- cpg_of_cluster %in% rownames(meth)
    grp      <- cluster_id_of[keep]
    cpg_ids  <- cpg_of_cluster[keep]

    if (length(grp) == 0) {
        mean_abs_cor_expr <- 0
    } else {
        sub_meth <- meth[cpg_ids, , drop = FALSE]
        avg_mat  <- rowsum(sub_meth, group = grp) / as.vector(table(grp))

        cluster_cors <- suppressWarnings(
            cor(t(avg_mat), expr, method = "spearman", use = "pairwise.complete.obs")[, 1]
        )
        mean_abs_cor_expr <- mean(abs(cluster_cors), na.rm = TRUE)
        if (!is.finite(mean_abs_cor_expr)) mean_abs_cor_expr <- 0
    }

    score <- modularity * (1 - (0.1 * fragmentation)) + j * mean_abs_cor_expr

    list(
        key      = key,
        score    = score,
        clusters = clusters,
        diag     = data.frame(
            K                 = key,
            n_clusters        = n_clusters,
            modularity        = modularity,
            fragmentation     = fragmentation,
            mean_abs_cor_expr = round(mean_abs_cor_expr, 3),
            score             = round(score, 4)
        )
    )
}

res_by_perc <- parallel::mclapply(
    mc.cores = C,
    X = percentiles_list,
    FUN = function(percentile_denom) {

        probs <- seq(0, 1, (1/percentile_denom))

        thresholds <- matrix(NA_real_, nrow = n_cpg, ncol = length(flex_points_list))
        for (i in seq_len(n_cpg)) {
            if (abs(row_medians[i]) < 0.1) next
            q <- quantile_from_sorted(row_sorted_asc[[i]], probs)
            d <- sort(diff(q), decreasing = TRUE)
            thresholds[i, ] <- d[flex_points_list]
        }

        out <- lapply(seq_along(flex_points_list), function(f) {
            KNN <- lapply(seq_len(n_cpg), function(i) {
                if (is.na(thresholds[i, f])) return(setNames(list(i), NULL)[[1]])
                thr <- thresholds[i, f]
                # count of values > thr via binary search on the ascending
                # sorted row, then slice that many indices off the
                # descending order -- replaces the O(n_cpg) which() scan
                n_above <- n_cpg - findInterval(thr, row_sorted_asc[[i]])
                neighbors.data <- if (n_above > 0) row_order_desc[[i]][seq_len(n_above)] else integer(0)
                ni <- intersect(neighbors.spatial.list[[i]], neighbors.data)
                if (length(ni) > 0) return(ni)
                NULL
            })
            names(KNN) <- coords$CpGID
            r <- connect_KNN_igraph(KNN)

            key <- paste0(percentile_denom, "_", flex_points_list[f])
            scored <- score_grid_point(r$clusters, r$graph, r$membership, key)

            rm(r, KNN)
            scored
        })

        rm(thresholds)
        out
    })

flat <- unlist(res_by_perc, recursive = FALSE)
rm(res_by_perc); invisible(gc(verbose = FALSE))

diag_df <- do.call(rbind, lapply(flat, `[[`, "diag"))

scores <- vapply(flat, function(x) {
    if (is.null(x$score) || !is.finite(x$score)) NA_real_ else x$score
}, numeric(1))

if (all(is.na(scores)))
  stop(sprintf("[%s] all %d grid points scored NA.", TxID, length(scores)))

if (all(is.na(scores))) {
    # clustering found no usable grouping -> each CpG is its own cluster
    cat(sprintf("[%s] no clustering (all grid points NA) - writing %d singleton CpGs\n",
                TxID, ncol(DC)))
    clusters_map <- data.frame(CpGID = colnames(DC), cluster_id = seq_len(ncol(DC)))
    write.table(clusters_map, file = output_path_map,
                row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
    pdf(output_path_plot); plot.new(); text(0.5, 0.5, "no clusters - singletons"); dev.off()
    quit(save = "no", status = 0)
}

best_i        <- which.max(scores)
best_k        <- flat[[best_i]]$key
KNN.connected <- flat[[best_i]]$clusters

cat(sprintf("[%s] best grid point: %s (%d clusters, score %.4f)\n",
            TxID, best_k, length(KNN.connected), scores[best_i]))

rm(flat, neighbors.spatial.list); invisible(gc(verbose = FALSE))

# ---- assign cluster ids ----
# CHANGE: we only ever needed rownames(DC) + a cluster label per CpG here --
# no need to carry a whole abs-value data.frame just to attach one column.
cluster_id_vec <- setNames(rep(NA_character_, n_cpg), rownames(DC))
for (idx in seq_along(KNN.connected)) {
    idxs <- KNN.connected[[idx]]
    cluster_id_vec[idxs] <- names(KNN.connected)[idx]
}

# ---- write the cluster map FIRST (the heatmap is the memory-hungry part) ----
clusters_map <- data.frame(CpGID = rownames(DC), cluster_id = unname(cluster_id_vec))
write.table(clusters_map, file = output_path_map, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
cat(sprintf("[%s] written: %s\n", TxID, output_path_map))

rm(DC); invisible(gc(verbose = FALSE))

# ---- heatmap ----
# CHANGE: re-read the (signed) matrix here instead of holding a
# dataset_plotting copy alongside DC for the whole script. This is the only
# point that needs signed values (for the blue/white/red color scale) --
# the cheap I/O here is a good trade for not carrying an extra n_cpg x n_cpg
# matrix through the entire grid search.
dataset_plotting <- read.delim(input_path, check.names = FALSE, header = TRUE, sep = "\t")
dataset_plotting$cluster_id <- cluster_id_vec[rownames(dataset_plotting)]

cluster_id <- dataset_plotting$cluster_id

mat <- as.matrix(dataset_plotting[, colnames(dataset_plotting) != "cluster_id"])
cpg_ids    <- rownames(mat)
start_pos <- as.numeric(coords$start)[match(cpg_ids, coords$CpGID)]
ord <- order(start_pos)

mat        <- mat[ord, ord]
cluster_id <- cluster_id[ord]
cluster_id <- as.character(cluster_id)

cluster_levels <- sort(as.numeric(unique(cluster_id)))
cluster_factor <- factor(cluster_id, levels = cluster_levels)

show_names <- nrow(mat) <= 300      # labels are unreadable and costly beyond this

ha_row <- rowAnnotation(
    cluster = cluster_factor,
    col = list(cluster = structure(rainbow(length(unique(cluster_factor))), names = levels(cluster_factor))),
    show_legend = FALSE
)

ha_col <- HeatmapAnnotation(
    cluster = cluster_factor,
    col = list(cluster = structure(rainbow(length(unique(cluster_factor))), names = levels(cluster_factor)))
)

col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

ht <- Heatmap(
    mat,
    name = "Correlation",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_split = cluster_factor,
    column_split = cluster_factor,
    top_annotation = ha_col,
    left_annotation = ha_row,
    show_row_names = show_names,
    show_column_names = show_names,
    column_names_gp = grid::gpar(fontsize = 4),
    row_names_gp = grid::gpar(fontsize = 4),
    col = col_fun,
    border = TRUE,
    row_title_rot = 90,
    column_title_rot = 90
)

pdf(output_path_plot, width = 10, height = 10)
draw(ht)
dev.off()

cat(sprintf("[%s] written: %s\n", TxID, output_path_plot))