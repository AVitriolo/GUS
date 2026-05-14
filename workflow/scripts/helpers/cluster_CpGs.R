cluster_CpGs <- function(params_grid, corr, coords, C){

    KNNs_list <- parallel::mclapply(
        mc.cores = C, 
        X = seq_len(nrow(params_grid)), 
        FUN = function(idx_out){

            percentile_denom <- params_grid$percentiles[idx_out]
            flex_point_idx <- params_grid$flex_points[idx_out]

            KNN <- parallel::mclapply(
                        mc.cores = C,
                        X = 1:nrow(coords), 
                        FUN = function(idx_in){

                            row <- unlist(corr[idx_in,])
                            percentiles <- quantile(row, probs = seq(0, 1, (1/percentile_denom)))

                            if (abs(median(row)) < 0.1){
                            neighbors.data <- idx_in
                            } else {
                            flex_point <- sort(diff(percentiles), decreasing = T)[flex_point_idx]
                            neighbors.data <- which(row >= flex_point)
                            }

                            neighbors.spatial <- find_neighbors(coords, idx_in, 1, T)
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
        data.frame(
            K             = k,
            n_edges       = igraph::ecount(r$graph),
            n_clusters    = length(unique(r$membership)),
            modularity    = round(igraph::modularity(r$graph, r$membership), 3),
            fragmentation = round(length(unique(r$membership)) / igraph::vcount(r$graph), 3)
        )
    }))

    diag_df$score <- diag_df$modularity * (1 - (0.1 * diag_df$fragmentation))

    best_k   <- diag_df$K[which.max(diag_df$score)]
    best_res <- KNNs_list[[best_k]]
    KNN.connected <- best_res$clusters

    corr[,"cluster_id"] <- NA

    for(idx in seq(KNN.connected)){
    idxs <- KNN.connected[[idx]]
    corr[idxs, "cluster_id"] <- names(KNN.connected)[idx]
    }

    n_clusters <- length(unique(corr$cluster_id))

    to.return <- list(ranking = diag_df, corr = corr, n_clusters = n_clusters)

    return(to.return)

}