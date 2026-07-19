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

# ---- load methylation matrix from rse for this TxID ----
rse  <- HDF5Array::loadHDF5SummarizedExperiment(input_path_rse)
meth <- SummarizedExperiment::assay(rse, "M")
meth <- as.matrix(meth)

# ---- LEAKAGE-SAFE: the j-term is scored on TRAIN samples only ----
# mean_abs_cor_expr selects the clustering parameters, so it must not see test
# methylation or test expression.
common_samples <- Reduce(intersect, list(colnames(meth), names(expr), TRAIN))
cat(sprintf("[%s] j-term scored on %d TRAIN samples\n", TxID, length(common_samples)))
if (j != 0 && length(common_samples) < 3)
  stop(sprintf("[%s] only %d train samples for the j-term.", TxID, length(common_samples)))

# =============================================================================
# Clustering grid search
#
# Same computation as before, reorganised to avoid redundant work:
#   - find_neighbors() depends only on (coords, i), so it is computed ONCE per
#     CpG instead of once per (CpG, percentile, flex_point).
#   - quantile()/diff() depend only on (CpG, percentile), so all flex_points
#     are read off the SAME sorted diff vector.
#   - the correlation matrix is used as a matrix, not a data.frame.
# Results are identical to the nested-loop version.
# =============================================================================
DC     <- as.matrix(dataset_clustering)
n_cpg  <- nrow(coords)

neighbors.spatial.list <- lapply(seq_len(n_cpg), function(i) find_neighbors(coords, i, 1, TRUE))
row_medians            <- matrixStats::rowMedians(DC)

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

        lapply(seq_along(flex_points_list), function(f) {
            KNN <- lapply(seq_len(n_cpg), function(i) {
                neighbors.data <- if (is.na(thresholds[i, f])) i else which(DC[i, ] >= thresholds[i, f])
                ni <- intersect(neighbors.spatial.list[[i]], neighbors.data)
                if (length(ni) > 0) return(ni)
                NULL
            })
            names(KNN) <- coords$CpGID
            connect_KNN_igraph(KNN)
        })
    })

# flatten to "<percentile>_<flex_point>" keys
KNNs_list <- list()
for (pi in seq_along(percentiles_list)) {
    for (fi in seq_along(flex_points_list)) {
        key <- paste0(percentiles_list[pi], "_", flex_points_list[fi])
        KNNs_list[[key]] <- res_by_perc[[pi]][[fi]]
    }
}
rm(res_by_perc); invisible(gc(verbose = FALSE))

# ---- score each grid point ----
diag_df <- do.call(rbind, lapply(names(KNNs_list), function(k) {
    r <- KNNs_list[[k]]

    n_edges       <- igraph::ecount(r$graph)
    n_clusters    <- length(unique(r$membership))
    modularity    <- round(igraph::modularity(r$graph, r$membership), 3)
    fragmentation <- round(n_clusters / igraph::vcount(r$graph), 3)

    clusters <- r$clusters
    # common_samples is TRAIN-only (defined above)

    cluster_cors <- sapply(clusters, function(idxs) {
        cpg_ids <- coords$CpGID[idxs]
        cpg_ids <- intersect(cpg_ids, rownames(meth))
        if (length(cpg_ids) == 0) return(NA)
        cluster_avg <- colMeans(meth[cpg_ids, common_samples, drop = FALSE], na.rm = TRUE)
        suppressWarnings(cor(cluster_avg, expr[common_samples], method = "spearman", use = "complete.obs"))
    })

    mean_abs_cor_expr <- mean(abs(cluster_cors), na.rm = TRUE)

    data.frame(
        K                 = k,
        n_edges           = n_edges,
        n_clusters        = n_clusters,
        modularity        = modularity,
        fragmentation     = fragmentation,
        mean_abs_cor_expr = round(mean_abs_cor_expr, 3)
    )
}))

diag_df$score <- diag_df$modularity * (1 - (0.1 * diag_df$fragmentation)) + j * diag_df$mean_abs_cor_expr

best_k   <- diag_df$K[which.max(diag_df$score)]
best_res <- KNNs_list[[best_k]]
KNN.connected <- best_res$clusters

cat(sprintf("[%s] best grid point: %s (%d clusters)\n",
            TxID, best_k, length(KNN.connected)))

# ---- assign cluster ids ----
dataset_clustering[,"cluster_id"] <- NA
dataset_plotting[,"cluster_id"]   <- NA

for(idx in seq(KNN.connected)){
    idxs <- KNN.connected[[idx]]
    dataset_clustering[idxs, "cluster_id"] <- names(KNN.connected)[idx]
    dataset_plotting[idxs, "cluster_id"]   <- names(KNN.connected)[idx]
}

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
    show_row_names = TRUE,
    show_column_names = TRUE,
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

clusters_map <- data.frame(CpGID = rownames(dataset_clustering), cluster_id = dataset_clustering$cluster_id)
write.table(clusters_map, file = output_path_map, row.names = F, col.names = F, quote = F, sep = "\t")

cat(sprintf("[%s] written: %s\n", TxID, output_path_map))