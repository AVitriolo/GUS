source("workflow/scripts/helpers/find_neighbors.R")
source("workflow/scripts/helpers/connect_KNN_igraph.R")

options(scipen=999)

library(FNN)
library(ComplexHeatmap)

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

input_path                     <- args$input_path
input_path_rse                 <- args$input_path_rse
input_path_counts              <- args$input_path_counts
C                               <- as.integer(args$C)
j                               <- as.numeric(args$j)
output_path_plot               <- args$output_path_plot
output_path_map                <- args$output_path_map

parts <- unlist(strsplit(basename(input_path), split = "_"))
TxID <- parts[length(parts)]

dataset <- read.delim(input_path, check.names = F, header = T, sep = "\t")

dataset_clustering <- abs(dataset)
dataset_plotting   <- dataset

percentiles_list <- seq(from = 10, to = 1000, length.out = 100)
flex_points_list <- seq(from = 1, to = 3, by = 1)

params_grid <- expand.grid(percentiles = percentiles_list, flex_points = flex_points_list)

CpGIDs <- colnames(dataset_clustering)
chrs   <- unlist(lapply(CpGIDs, function(x) strsplit(x, "_")[[1]][1]))
starts <- unlist(lapply(CpGIDs, function(x) strsplit(x, "_")[[1]][2]))
coords <- data.frame(chr = chrs, start = starts, end = starts, CpGID = CpGIDs)

# ---- load expression for this TxID ----
counts <- read.table(input_path_counts, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
expr   <- as.numeric(counts[TxID, ])
names(expr) <- colnames(counts)

# ---- load methylation matrix from rse for this TxID ----
rse  <- HDF5Array::loadHDF5SummarizedExperiment(input_path_rse)
meth <- SummarizedExperiment::assay(rse, "M")
meth <- as.matrix(meth)

KNNs_list <- parallel::mclapply(
    mc.cores = C,
    X = seq_len(nrow(params_grid)),
    FUN = function(idx_out){

        percentile_denom <- params_grid$percentiles[idx_out]
        flex_point_idx    <- params_grid$flex_points[idx_out]

        KNN <- parallel::mclapply(
            mc.cores = C,
            X = 1:nrow(coords),
            FUN = function(idx_in){

                row <- unlist(dataset_clustering[idx_in,])
                percentiles <- quantile(row, probs = seq(0, 1, (1/percentile_denom)))

                if (abs(median(row)) < 0.1){
                    neighbors.data <- idx_in
                } else {
                    flex_point <- sort(diff(percentiles), decreasing = T)[flex_point_idx]
                    neighbors.data <- which(row >= flex_point)
                }

                neighbors.spatial   <- find_neighbors(coords, idx_in, 1, T)
                neighbors.intersect <- intersect(neighbors.spatial, neighbors.data)

                if (length(neighbors.intersect) > 0){return(neighbors.intersect)}; return(NULL)
            })

        names(KNN) <- coords$CpGID
        KNN.connected <- connect_KNN_igraph(KNN)

        return(KNN.connected)
    })

names(KNNs_list) <- paste0(params_grid$percentiles, "_", params_grid$flex_points)

diag_df <- do.call(rbind, lapply(names(KNNs_list), function(k) {
    r <- KNNs_list[[k]]

    n_edges       <- igraph::ecount(r$graph)
    n_clusters    <- length(unique(r$membership))
    modularity    <- round(igraph::modularity(r$graph, r$membership), 3)
    fragmentation <- round(n_clusters / igraph::vcount(r$graph), 3)

    # mean abs correlation between cluster-averaged methylation and expression
    clusters <- r$clusters
    common_samples <- intersect(colnames(meth), names(expr))

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

# ---- score: unsupervised term + j * expression-correlation term ----
diag_df$score <- diag_df$modularity * (1 - (0.1 * diag_df$fragmentation)) + j * diag_df$mean_abs_cor_expr

best_k   <- diag_df$K[which.max(diag_df$score)]
best_res <- KNNs_list[[best_k]]
KNN.connected <- best_res$clusters

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