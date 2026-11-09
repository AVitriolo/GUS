# =============================================================================
# most_similar_n_regions_min_dist()
# As most_similar_n_regions.R (PCA of the centred, scaled feature matrix,
# first n_pcs PCs, kNN), but a candidate on the same chromosome closer than
# min_dist bp (gap between the two windows; overlapping and identical windows
# included) is skipped and the next most similar region is taken, so a region
# is not compared with its own locus.
#   x          regions x features matrix (rownames = region names)
#   chr/start/end  coordinates of the rows of x
# Returns a named list (rownames of x) of neighbour row indices, most similar first.
# =============================================================================
most_similar_n_regions_min_dist <- function(x, chr, start, end, n = 100, n_pcs = 20,
                                            min_dist = 1000, num_cores = 8) {
    mu <- colMeans(x); sdv <- sqrt(pmax(colMeans(x^2) - mu^2, 0))
    if (any(sdv == 0)) message(sprintf("dropping %d constant features", sum(sdv == 0)))
    x <- x[, sdv > 0, drop = FALSE]
    pca <- prcomp(x, center = TRUE, scale. = TRUE, rank. = n_pcs)$x
    rm(x); invisible(gc())
    message("donePCA")

    N <- nrow(pca); chr_i <- as.integer(factor(chr)); st <- start; en <- end
    nb   <- matrix(NA_integer_, N, n)
    todo <- seq_len(N); k <- min(N, 3L * n); n_skipped <- numeric(0)
    while (length(todo)) {
        chunks <- split(todo, ceiling(seq_along(todo) / min(20000, max(1, floor(1e8 / k)))))   # ~1 GB per worker
        R <- do.call(rbind, parallel::mclapply(chunks, function(ch)
            RANN::nn2(pca, pca[ch, , drop = FALSE], k = k)$nn.idx, mc.cores = num_cores))
        gap  <- pmax(matrix(st[R], nrow(R)) - en[todo], st[todo] - matrix(en[R], nrow(R)))
        near <- R == todo | (matrix(chr_i[R], nrow(R)) == chr_i[todo] & gap < min_dist)
        if (!length(n_skipped)) n_skipped <- rowSums((near & R != todo)[, seq_len(min(k, n + 1L)), drop = FALSE])
        R[near] <- NA_integer_
        ok <- rowSums(!is.na(R)) >= n | k >= N
        if (any(ok))
            nb[todo[ok], ] <- t(apply(R[ok, , drop = FALSE], 1, function(v) c(v[!is.na(v)], rep(NA_integer_, n))[seq_len(n)]))
        message(sprintf("k = %d: %d regions done, %d need more candidates", k, sum(ok), sum(!ok)))
        todo <- todo[!ok]; k <- min(N, 4L * k)
        rm(R, gap, near); invisible(gc())
    }
    message(sprintf("neighbours within %d bp skipped: in the first %d candidates, mean %.2f, >0 for %.1f%% of regions",
                    min_dist, n, mean(n_skipped), 100 * mean(n_skipped > 0)))
    message("doneknn")
    setNames(lapply(seq_len(N), function(i) { v <- nb[i, ]; v[!is.na(v)] }), rownames(pca))
}
