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

dataset <- read.delim(input_path, check.names = F, header = T, sep = "\t")

dataset_clustering <- abs(dataset)
dataset_plotting   <- dataset

percentiles_list <- seq(from = 10, to = 1000, length.out = 100)
flex_points_list <- seq(from = 1, to = 3, by = 1)

CpGIDs <- colnames(dataset_clustering)
chrs   <- unlist(lapply(CpGIDs, function(x) strsplit(x, "_")[[1]][1]))
starts <- unlist(lapply(CpGIDs, function(x) strsplit(x, "_")[[1]][2]))
coords <- data.frame(chr = chrs, start = starts, end = starts, CpGID = CpGIDs)

# ---- load expression for this TxID (TRAIN counts) ----
counts <- read.table(input_path_counts, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
expr   <- as.numeric(counts[TxID, ])
names(expr) <- colnames(counts)
rm(counts); invisible(gc(verbose = FALSE))

# ---- load methylation matrix from rse for this TxID ----
rse  <- HDF5Array::loadHDF5SummarizedExperiment(input_path_rse)
meth <- SummarizedExperiment::assay(rse, "M")
meth <- as.matrix(meth)
rm(rse); invisible(gc(verbose = FALSE))

# ---- LEAKAGE-SAFE: the j-term is scored on TRAIN samples only ----
# mean_abs_cor_expr selects the clustering parameters, so it must not see test
# methylation or test expression.
common_samples <- Reduce(intersect, list(colnames(meth), names(expr), TRAIN))
cat(sprintf("[%s] j-term scored on %d TRAIN samples\n", TxID, length(common_samples)))
if (j != 0 && length(common_samples) < 3)
  stop(sprintf("[%s] only %d train samples for the j-term.", TxID, length(common_samples)))

# keep only what the j-term needs, so the big matrix isn't carried into the workers
meth <- meth[, common_samples, drop = FALSE]
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
#     so only the cluster lists survive (the previous version held all 300
#     graphs at once, which is what triggered the OOM kill)
# =============================================================================
DC     <- as.matrix(dataset_clustering)
n_cpg  <- nrow(coords)

cat(sprintf("[%s] %d CpGs, %d grid points\n",
            TxID, n_cpg, length(percentiles_list) * length(flex_points_list)))

neighbors.spatial.list <- lapply(seq_len(n_cpg), function(i) find_neighbors(coords, i, 1, TRUE))
row_medians            <- matrixStats::rowMedians(DC)

score_grid_point <- function(clusters, graph, membership, key) {
    n_clusters    <- length(unique(membership))
    modularity    <- igraph::modularity(graph, membership)
    if (!is.finite(modularity)) modularity <- 0      # no edges -> undefined
    modularity    <- round(modularity, 3)
    fragmentation <- round(n_clusters / igraph::vcount(graph), 3)

    cluster_cors <- sapply(clusters, function(idxs) {
        cpg_ids <- intersect(coords$CpGID[idxs], rownames(meth))
        if (length(cpg_ids) == 0) return(NA)
        cluster_avg <- colMeans(meth[cpg_ids, , drop = FALSE], na.rm = TRUE)
        suppressWarnings(cor(cluster_avg, expr, method = "spearman", use = "complete.obs"))
    })
    mean_abs_cor_expr <- mean(abs(cluster_cors), na.rm = TRUE)
    if (!is.finite(mean_abs_cor_expr)) mean_abs_cor_expr <- 0

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

        # thresholds[i, f] : flex threshold for CpG i at flex_point f
        # NA marks the "singleton" case (|median| < 0.1)
        thresholds <- matrix(NA_real_, nrow = n_cpg, ncol = length(flex_points_list))
        for (i in seq_len(n_cpg)) {
            if (abs(row_medians[i]) < 0.1) next
            d <- sort(diff(quantile(DC[i, ], probs = probs)), decreasing = TRUE)
            thresholds[i, ] <- d[flex_points_list]
        }

        out <- lapply(seq_along(flex_points_list), function(f) {
            KNN <- lapply(seq_len(n_cpg), function(i) {
                neighbors.data <- if (is.na(thresholds[i, f])) i else which(DC[i, ] >= thresholds[i, f])
                ni <- intersect(neighbors.spatial.list[[i]], neighbors.data)
                if (length(ni) > 0) return(ni)
                NULL
            })
            names(KNN) <- coords$CpGID
            r <- connect_KNN_igraph(KNN)

            key <- paste0(percentile_denom, "_", flex_points_list[f])
            scored <- score_grid_point(r$clusters, r$graph, r$membership, key)

            rm(r, KNN)                       # drop the igraph object immediately
            scored
        })

        rm(thresholds); invisible(gc(verbose = FALSE))
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

best_i        <- which.max(scores)
best_k        <- flat[[best_i]]$key
KNN.connected <- flat[[best_i]]$clusters


cat(sprintf("[%s] best grid point: %s (%d clusters, score %.4f)\n",
            TxID, best_k, length(KNN.connected), scores[best_i]))

rm(flat, DC, neighbors.spatial.list); invisible(gc(verbose = FALSE))

# ---- assign cluster ids ----
dataset_clustering[,"cluster_id"] <- NA
dataset_plotting[,"cluster_id"]   <- NA

for(idx in seq(KNN.connected)){
    idxs <- KNN.connected[[idx]]
    dataset_clustering[idxs, "cluster_id"] <- names(KNN.connected)[idx]
    dataset_plotting[idxs, "cluster_id"]   <- names(KNN.connected)[idx]
}

# ---- write the cluster map FIRST (the heatmap is the memory-hungry part) ----
clusters_map <- data.frame(CpGID = rownames(dataset_clustering), cluster_id = dataset_clustering$cluster_id)
write.table(clusters_map, file = output_path_map, row.names = F, col.names = F, quote = F, sep = "\t")
cat(sprintf("[%s] written: %s\n", TxID, output_path_map))

# ---- heatmap ----
cluster_id <- dataset_plotting$cluster_id

mat <- as.matrix(dataset_plotting[, colnames(dataset_plotting) != "cluster_id"])
cpg_ids    <- rownames(mat)
start_pos  <- as.numeric(sapply(strsplit(cpg_ids, "_"), function(x) x[2]))
ord        <- order(start_pos)

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