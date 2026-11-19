# Fast, result-identical replacement for
#   compute_pairwise_correlation_split.R  +  merge_CpGs_split_new.R
# (the original scripts are untouched; this one reads the RSE directly).
#
# What the original clustering actually uses from the CpG x CpG |Spearman| matrix DC:
#   * per CpG i: the median of row i, and the 300 thresholds
#       sort(diff(quantile(DC[i, ], seq(0, 1, 1/denom))), decreasing = TRUE)[1:3]
#     for denom in seq(10, 1000, length.out = 100)
#   * DC[i, i-1], DC[i, i], DC[i, i+1]   (find_neighbors(K = 1, return_self = TRUE))
# so the matrix is computed in row blocks (ranks once + BLAS crossprod), these numbers
# are extracted, and the block is discarded: no multi-GB text file, no full matrix in RAM.
#
# Exactness. The original values are pairwise Spearman (C_cor on per-pair ranks),
# written by write.table (15 significant digits) and read back. The BLAS values differ
# from those by < 1e-15 (measured: <= 7e-16). Every clustering decision is a comparison
#   |median_i| < 0.1        and        DC[i, nb] >= threshold_i,g
# and perturbing every DC value by <= d moves medians/quantiles by <= d and the
# thresholds (differences of quantiles) by <= 2d. So when every comparison of row i has
# a margin > TOL (1e-12 >> 3d) the decisions -- hence KNN, clusters, scores, best grid
# point and the map -- are identical. Rows with any comparison inside TOL are recomputed
# EXACTLY as the original does (per-pair C_cor, write.table-equivalent 15-digit rounding,
# original quantile code) before the grid search.
# The grid search reproduces find_neighbors + connect_KNN_igraph exactly (KNN sets are
# subsets of {i-1, i, i+1}, so set overlaps reduce to the edge rules below) and scores
# with the original score_grid_point (copied verbatim).
#
# Heatmap: identical layout; for genes with > PLOT_MAX CpGs it shows an evenly spaced
# subsample of PLOT_MAX CpGs (the original drew every CpG). Only the figure differs.

options(scipen = 999)
library(ComplexHeatmap)

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
input_path_rse     <- args$input_path_rse
input_path_counts  <- args$input_path_counts
train_samples_path <- args$train_samples
C                  <- as.integer(args$C)
j                  <- as.numeric(args$j)
output_path_plot   <- args$output_path_plot
output_path_map    <- args$output_path_map
output_path_diag   <- args$output_path_diag            # optional: grid diagnostics (for validation)
TxID               <- args$tx_id
BLOCK    <- 256L
TOL      <- if (is.null(args$tol)) 1e-12 else as.numeric(args$tol)   # --tol only to test the exact path
PLOT_MAX <- 2000L

TRAIN <- readLines(train_samples_path)
dir.create(dirname(output_path_map),  recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_path_plot), recursive = TRUE, showWarnings = FALSE)
t_start <- Sys.time()

# =============================================================================
# Inputs, exactly as the two original scripts build them
# =============================================================================
rse  <- HDF5Array::loadHDF5SummarizedExperiment(input_path_rse)
meth_full <- as.matrix(SummarizedExperiment::assay(rse, "M"))       # CpGs x samples

# --- correlation input (compute_pairwise_correlation_split.R) ---
coord_order <- GenomicRanges::mcols(GenomicRanges::sort(SummarizedExperiment::rowRanges(rse)))$CpGID
X <- t(meth_full)
X <- X[rownames(X) %in% TRAIN, , drop = FALSE]                       # train samples, RSE order
cat(sprintf("[%s] Samples: %d -> %d train\n", TxID, ncol(meth_full), nrow(X)))
if (nrow(X) < 3) stop(sprintf("Only %d train samples - cannot compute correlations.", nrow(X)))
coord_order <- coord_order[coord_order %in% colnames(X)]
X <- X[, coord_order, drop = FALSE]
if (anyNA(X))
  stop(sprintf("[%s] NA in train methylation: the fast path assumes complete data; use the original rules.", TxID))

CpGIDs <- coord_order
n_cpg  <- length(CpGIDs)
parts_split <- strsplit(CpGIDs, "_", fixed = TRUE)
coords <- data.frame(chr   = vapply(parts_split, `[`, character(1), 1),
                     start = vapply(parts_split, `[`, character(1), 2))
coords$end <- coords$start; coords$CpGID <- CpGIDs

# --- j-term input (merge_CpGs_split_new.R) ---
counts <- read.table(input_path_counts, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
expr   <- as.numeric(counts[TxID, ]); names(expr) <- colnames(counts)
rm(counts)
common_samples <- Reduce(intersect, list(colnames(rse), names(expr), TRAIN))
cat(sprintf("[%s] j-term scored on %d TRAIN samples\n", TxID, length(common_samples)))
if (j != 0 && length(common_samples) < 3)
  stop(sprintf("[%s] only %d train samples for the j-term.", TxID, length(common_samples)))
meth <- meth_full[, common_samples, drop = FALSE]
expr <- expr[common_samples]
rm(rse, meth_full); invisible(gc(verbose = FALSE))

# =============================================================================
# Correlations in row blocks -> per-row medians, thresholds, neighbour values
# =============================================================================
percentiles_list <- seq(from = 10, to = 1000, length.out = 100)
flex_points_list <- seq(from = 1, to = 3, by = 1)
n_grid <- length(percentiles_list) * length(flex_points_list)
cat(sprintf("[%s] %d CpGs, %d grid points\n", TxID, n_cpg, n_grid))

# quantile index tables, built with the original expressions
quantile_from_sorted <- function(sorted_asc, probs) {
    n  <- length(sorted_asc)
    h  <- (n - 1) * probs + 1
    lo <- pmin(pmax(floor(h), 1), n)
    hi <- pmin(pmax(ceiling(h), 1), n)
    sorted_asc[lo] + (h - lo) * (sorted_asc[hi] - sorted_asc[lo])
}
qtab <- lapply(percentiles_list, function(d) {
    probs <- seq(0, 1, (1 / d))
    h  <- (n_cpg - 1) * probs + 1
    lo <- pmin(pmax(floor(h), 1), n_cpg)
    hi <- pmin(pmax(ceiling(h), 1), n_cpg)
    list(lo = lo, hi = hi, fr = h - lo)
})
LO <- unlist(lapply(qtab, `[[`, "lo")); HI <- unlist(lapply(qtab, `[[`, "hi"))
FR <- unlist(lapply(qtab, `[[`, "fr"))
qlen <- vapply(qtab, function(z) length(z$lo), 1L)
qend <- cumsum(qlen); qbeg <- qend - qlen + 1L
# consecutive pairs inside each percentile group (diff(q))
D2 <- unlist(lapply(seq_along(qlen), function(g) (qbeg[g] + 1L):qend[g]))
D1 <- D2 - 1L
dgrp <- rep(seq_along(qlen), qlen - 1L)
dcols <- split(seq_along(D2), dgrp)

# ranks once per CpG (cor() ranks each column with rank(); no NA -> identical)
Rk <- apply(X, 2L, rank)
if (!is.matrix(Rk)) Rk <- matrix(Rk, ncol = n_cpg)
Zc <- sweep(Rk, 2L, colMeans(Rk))
nrm <- sqrt(colSums(Zc^2))
const <- nrm == 0                                   # zero-variance CpGs -> NA in the original
has_na <- any(const) && n_cpg >= 2                  # original: NA -> 0 and diag -> 1 only if any NA
Z <- sweep(Zc, 2L, ifelse(const, 1, nrm), "/"); Z[, const] <- 0
rm(Zc)

has_prev <- c(FALSE, coords$chr[-1] == coords$chr[-n_cpg])
has_next <- c(coords$chr[-1] == coords$chr[-n_cpg], FALSE)
if (n_cpg == 1) { has_prev <- FALSE; has_next <- FALSE }

plot_idx <- if (n_cpg <= PLOT_MAX) seq_len(n_cpg) else unique(round(seq(1, n_cpg, length.out = PLOT_MAX)))

# corrected signed correlations of rows `rows` (approximate, BLAS)
signed_block <- function(rows) {
    if (n_cpg == 1) return(matrix(1, 1, 1))
    B <- crossprod(Z[, rows, drop = FALSE], Z)
    B[B > 1] <- 1; B[B < -1] <- -1
    if (has_na) {
        B[const[rows], ] <- 0; B[, const] <- 0
        B[cbind(seq_along(rows), rows)] <- 1
    }
    B
}

thr_from_sorted <- function(S) {      # S: rows x n sorted |r| -> rows x 300 thresholds
    Q <- S[, LO, drop = FALSE] + sweep(S[, HI, drop = FALSE] - S[, LO, drop = FALSE], 2L, FR, "*")
    Dm <- Q[, D2, drop = FALSE] - Q[, D1, drop = FALSE]
    out <- matrix(NA_real_, nrow(S), n_grid)
    rr <- seq_len(nrow(S))
    for (g in seq_along(dcols)) {
        sub <- Dm[, dcols[[g]], drop = FALSE]
        # sort(r, decreasing = TRUE)[1:3] per row, as 3 x (row max, drop one occurrence):
        # same values, duplicates included (diff() of sorted quantiles has no NA)
        for (f in flex_points_list) {
            out[, (g - 1) * 3 + f] <- matrixStats::rowMaxs(sub)
            if (f < max(flex_points_list)) sub[cbind(rr, max.col(sub, ties.method = "first"))] <- -Inf
        }
    }
    out
}

process_block <- function(rows) {
    A  <- signed_block(rows)
    sg <- if (length(plot_idx) && any(rows %in% plot_idx)) A[rows %in% plot_idx, plot_idx, drop = FALSE] else NULL
    A  <- abs(A)
    S  <- t(apply(A, 1L, sort)); if (n_cpg == 1) S <- matrix(A, 1, 1)
    med <- matrixStats::rowMedians(S)
    thr <- thr_from_sorted(S)
    thr[abs(med) < 0.1, ] <- NA
    k  <- seq_along(rows)
    vself <- A[cbind(k, rows)]
    vprev <- ifelse(has_prev[rows], A[cbind(k, pmax(rows - 1L, 1L))], NA)
    vnext <- ifelse(has_next[rows], A[cbind(k, pmin(rows + 1L, n_cpg))], NA)
    # smallest distance of any decision of the row from its cut-off
    marg <- abs(med - 0.1)
    for (v in list(vprev, vself, vnext)) {
        dd <- abs(thr - v)                     # NA where no threshold / no neighbour
        m  <- suppressWarnings(apply(dd, 1L, min, na.rm = TRUE))
        marg <- pmin(marg, m)
    }
    list(rows = rows, thr = thr, vprev = vprev, vself = vself, vnext = vnext,
         marg = marg, sg = sg)
}

blocks <- split(seq_len(n_cpg), ceiling(seq_len(n_cpg) / BLOCK))
t0 <- Sys.time()
res <- parallel::mclapply(blocks, process_block, mc.cores = C)
if (any(vapply(res, inherits, TRUE, "try-error"))) stop("block computation failed")
rows_all <- unname(unlist(lapply(res, `[[`, "rows")))
stopifnot(identical(rows_all, seq_len(n_cpg)))
thresholds <- do.call(rbind, lapply(res, `[[`, "thr"))
vprev <- unlist(lapply(res, `[[`, "vprev")); vself <- unlist(lapply(res, `[[`, "vself"))
vnext <- unlist(lapply(res, `[[`, "vnext")); marg <- unlist(lapply(res, `[[`, "marg"))
sg_rows <- do.call(rbind, lapply(res, `[[`, "sg"))
rm(res); invisible(gc(verbose = FALSE))
cat(sprintf("[%s] block pass: %.1f s\n", TxID, as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# ---- exact recomputation of rows with a decision inside TOL ----
exact_row <- function(i) {
    if (n_cpg == 1) { r <- 1 } else {
        r <- vapply(seq_len(n_cpg), function(k)
            .Call(stats:::C_cor, Rk[, max(i, k)], Rk[, min(i, k)], 1L, FALSE), numeric(1))
        if (has_na) { r[is.na(r)] <- 0; r[i] <- 1 }
    }
    # write.table(15 significant digits) + read.delim round trip
    r <- as.numeric(vapply(r, function(x) format(x, digits = 15), ""))
    abs(r)
}
redo <- which(marg <= TOL)
cat(sprintf("[%s] %d row(s) with a decision within %g of its cut-off -> exact recomputation\n",
            TxID, length(redo), TOL))
for (i in redo) {
    DCi <- exact_row(i)
    s   <- sort(DCi)
    th  <- rep(NA_real_, n_grid)
    if (!(abs(matrixStats::rowMedians(matrix(DCi, 1))) < 0.1)) {
        for (g in seq_along(percentiles_list)) {
            q <- quantile_from_sorted(s, seq(0, 1, (1 / percentiles_list[g])))
            d <- sort(diff(q), decreasing = TRUE)
            th[(g - 1) * 3 + 1:3] <- d[flex_points_list]
        }
    }
    thresholds[i, ] <- th
    vself[i] <- DCi[i]
    vprev[i] <- if (has_prev[i]) DCi[i - 1] else NA
    vnext[i] <- if (has_next[i]) DCi[i + 1] else NA
}
rm(Z); invisible(gc(verbose = FALSE))

# =============================================================================
# Grid search: find_neighbors + connect_KNN_igraph, exact, vectorised
# =============================================================================
score_grid_point <- function(clusters, graph, membership, key) {
    n_clusters    <- length(unique(membership))
    modularity    <- igraph::modularity(graph, membership)
    if (!is.finite(modularity)) modularity <- 0      # no edges -> undefined
    modularity    <- round(modularity, 3)
    fragmentation <- round(n_clusters / igraph::vcount(graph), 3)

    # ---- vectorized cluster-mean correlation, built off `clusters` (not `membership`) ----
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

# Same score, faster j-term. cor(t(avg_mat), expr, "spearman", "pairwise.complete.obs") loops
# in R over clusters and, per cluster, does complete.cases, rank(cluster), rank(expr) and
#   .Call(C_cor, rank(x), rank(y), 1L, FALSE)
# With no NA (all cases complete) this makes exactly that C call on exactly the same
# ranks (rowRanks "average" == rank()), so every correlation is bit-identical.
no_na_score <- !anyNA(meth) && !anyNA(expr)
rank_expr   <- rank(as.numeric(expr))
C_cor <- get("C_cor", envir = asNamespace("stats"))
score_grid_point_fast <- function(clusters, graph, membership, key) {
    if (!no_na_score) return(score_grid_point(clusters, graph, membership, key))
    n_clusters    <- length(unique(membership))
    modularity    <- igraph::modularity(graph, membership)
    if (!is.finite(modularity)) modularity <- 0
    modularity    <- round(modularity, 3)
    fragmentation <- round(n_clusters / igraph::vcount(graph), 3)

    cluster_id_of  <- rep(seq_along(clusters), times = lengths(clusters))
    cpg_of_cluster <- coords$CpGID[unlist(clusters)]
    keep    <- cpg_of_cluster %in% rownames(meth)
    grp     <- cluster_id_of[keep]
    cpg_ids <- cpg_of_cluster[keep]

    if (length(grp) == 0) {
        mean_abs_cor_expr <- 0
    } else {
        sub_meth <- meth[cpg_ids, , drop = FALSE]
        avg_mat  <- rowsum(sub_meth, group = grp) / as.vector(table(grp))
        rk <- matrixStats::rowRanks(avg_mat, ties.method = "average")
        cluster_cors <- suppressWarnings(vapply(seq_len(nrow(rk)), function(c)
            .Call(C_cor, rk[c, ], rank_expr, 1L, FALSE), numeric(1)))
        mean_abs_cor_expr <- mean(abs(cluster_cors), na.rm = TRUE)
        if (!is.finite(mean_abs_cor_expr)) mean_abs_cor_expr <- 0
    }
    score <- modularity * (1 - (0.1 * fragmentation)) + j * mean_abs_cor_expr
    list(key = key, score = score, clusters = clusters,
         diag = data.frame(K = key, n_clusters = n_clusters, modularity = modularity,
                           fragmentation = fragmentation,
                           mean_abs_cor_expr = round(mean_abs_cor_expr, 3), score = round(score, 4)))
}
VERIFY <- !is.null(args$verify_scores) && as.integer(args$verify_scores) == 1   # validation only

n  <- n_cpg
vp <- ifelse(is.na(vprev), -Inf, vprev); vn <- ifelse(is.na(vnext), -Inf, vnext)

knn_grid_point <- function(g) {
    t    <- thresholds[, g]
    isna <- is.na(t)                                  # |median| < 0.1 -> KNN = {i}
    t2   <- ifelse(isna, Inf, t)
    inP  <- !isna & has_prev & (vp >= t2)             # i-1 in KNN[[i]]
    inN  <- !isna & has_next & (vn >= t2)             # i+1 in KNN[[i]]
    inS  <- isna | (vself >= t2)                      # i   in KNN[[i]]
    inS[!inP & !inS & !inN] <- TRUE                   # NULL -> i (connect_KNN_igraph)
    # KNN[[i]] and KNN[[j]] (i < j) share an element only if
    #   j = i+1: {i} (inS[i] & inP[j}) or {i+1} (inN[i] & inS[j]);   j = i+2: {i+1} (inN[i] & inP[j])
    e1 <- if (n >= 2) which((inS[-n] & inP[-1]) | (inN[-n] & inS[-1])) else integer(0)
    e2 <- if (n >= 3) which(inN[1:(n - 2)] & inP[3:n]) else integer(0)
    ei <- c(e1, e2); ej <- c(e1 + 1L, e2 + 2L)
    o  <- order(ei, ej)
    edges <- as.integer(rbind(ei[o], ej[o]))
    gr <- igraph::graph(edges = edges, directed = FALSE, n = n)
    comps <- igraph::components(gr)
    memb  <- comps$membership
    # cluster = sorted unique union of the members' KNN sets
    w_s <- which(inS); w_p <- which(inP); w_n <- which(inN)
    el_m <- memb[c(w_s, w_p, w_n)]
    el_v <- c(w_s, w_p - 1L, w_n + 1L)
    keep <- !duplicated(cbind(el_m, el_v))
    el_m <- el_m[keep]; el_v <- el_v[keep]
    o <- order(el_m, el_v)
    clusters <- split(as.integer(el_v[o]), factor(el_m[o], levels = sort(unique(memb))))
    list(clusters = clusters, graph = gr, membership = memb, edges = edges)
}

grid_label <- function(g) {
    p <- (g - 1) %/% length(flex_points_list) + 1; f <- (g - 1) %% length(flex_points_list) + 1
    paste0(percentiles_list[p], "_", flex_points_list[f])
}

# Pass 1: fingerprint each grid point's graph (edge list) + clusters. Grid points with the
# same fingerprint have the same graph, membership and clusters, hence the same score.
t0 <- Sys.time()
fp <- unlist(parallel::mclapply(seq_len(n_grid), function(g) {
    r <- knn_grid_point(g)
    digest::digest(list(r$edges, unname(r$clusters), names(r$clusters)), algo = "sha1")
}, mc.cores = C))
rep_of <- match(fp, fp)                       # first grid point with the same fingerprint
reps   <- unique(rep_of)
cat(sprintf("[%s] %d distinct clusterings among %d grid points\n", TxID, length(reps), n_grid))

# Pass 2: score each distinct clustering once
scored_reps <- parallel::mclapply(reps, function(g) {
    r <- knn_grid_point(g)
    s <- score_grid_point_fast(r$clusters, r$graph, r$membership, grid_label(g))
    if (VERIFY) {
        s0 <- score_grid_point(r$clusters, r$graph, r$membership, grid_label(g))
        if (!identical(s0$score, s$score) || !identical(s0$diag, s$diag))
            stop(sprintf("score mismatch at grid point %s: %.17g vs %.17g", grid_label(g), s0$score, s$score))
    }
    s
}, mc.cores = C)
if (any(vapply(scored_reps, inherits, TRUE, "try-error"))) stop(as.character(scored_reps[vapply(scored_reps, inherits, TRUE, "try-error")][[1]]))
names(scored_reps) <- as.character(reps)
flat <- lapply(seq_len(n_grid), function(g) {
    s <- scored_reps[[as.character(rep_of[g])]]
    s$key <- grid_label(g); s$diag$K <- s$key
    s
})
rm(scored_reps); invisible(gc(verbose = FALSE))
if (VERIFY) cat(sprintf("[%s] verify_scores: fast scores bit-identical to score_grid_point\n", TxID))
cat(sprintf("[%s] grid search: %.1f s\n", TxID, as.numeric(difftime(Sys.time(), t0, units = "secs"))))

diag_df <- do.call(rbind, lapply(flat, `[[`, "diag"))
if (!is.null(output_path_diag))
    write.table(diag_df, output_path_diag, sep = "\t", quote = FALSE, row.names = FALSE)

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
rm(flat); invisible(gc(verbose = FALSE))

# ---- assign cluster ids (original code) ----
cluster_id_vec <- setNames(rep(NA_character_, n_cpg), CpGIDs)
for (idx in seq_along(KNN.connected)) {
    idxs <- KNN.connected[[idx]]
    cluster_id_vec[idxs] <- names(KNN.connected)[idx]
}
clusters_map <- data.frame(CpGID = CpGIDs, cluster_id = unname(cluster_id_vec))
write.table(clusters_map, file = output_path_map, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
cat(sprintf("[%s] written: %s\n", TxID, output_path_map))

# =============================================================================
# Heatmap (original layout; evenly spaced subsample above PLOT_MAX CpGs)
# =============================================================================
mat <- sg_rows; dimnames(mat) <- list(CpGIDs[plot_idx], CpGIDs[plot_idx])
cluster_id <- as.character(cluster_id_vec[plot_idx])     # plot_idx is already in coordinate order
cluster_levels <- sort(as.numeric(unique(cluster_id)))
cluster_factor <- factor(cluster_id, levels = cluster_levels)
show_names <- nrow(mat) <= 300
col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
MAX_CLUSTERS_FOR_SPLIT <- 500
n_final_clusters <- length(unique(cluster_id))
cat(sprintf("[%s] %d clusters over %d plotted CpGs (of %d) -- building heatmap...\n",
            TxID, n_final_clusters, nrow(mat), n_cpg))

PLOT <- is.null(args$plot) || as.integer(args$plot) == 1
if (!PLOT) {
    pdf(output_path_plot); plot.new()
    text(0.5, 0.5, sprintf("heatmap disabled (clustering_heatmap: 0)\n%d CpGs, %d clusters", n_cpg, n_final_clusters)); dev.off()
} else if (nrow(mat) < 2) {
    pdf(output_path_plot); plot.new()
    text(0.5, 0.5, sprintf("%d CpG(s) - no heatmap", nrow(mat))); dev.off()
} else {
    ttl <- if (n_cpg > PLOT_MAX) sprintf("%d of %d CpGs (evenly spaced)", nrow(mat), n_cpg) else NULL
    if (n_final_clusters > MAX_CLUSTERS_FOR_SPLIT) {
        ht <- Heatmap(mat, name = "Correlation", cluster_rows = FALSE, cluster_columns = FALSE,
                      show_row_names = FALSE, show_column_names = FALSE, col = col_fun, border = TRUE)
    } else {
        cols <- structure(rainbow(length(unique(cluster_factor))), names = levels(cluster_factor))
        ha_row <- rowAnnotation(cluster = cluster_factor, col = list(cluster = cols), show_legend = FALSE)
        ha_col <- HeatmapAnnotation(cluster = cluster_factor, col = list(cluster = cols))
        ht <- Heatmap(mat, name = "Correlation", cluster_rows = FALSE, cluster_columns = FALSE,
                      row_split = cluster_factor, column_split = cluster_factor,
                      top_annotation = ha_col, left_annotation = ha_row,
                      show_row_names = show_names, show_column_names = show_names,
                      column_names_gp = grid::gpar(fontsize = 4), row_names_gp = grid::gpar(fontsize = 4),
                      col = col_fun, border = TRUE, row_title_rot = 90, column_title_rot = 90)
    }
    pdf(output_path_plot, width = 10, height = 10)
    draw(ht, column_title = ttl)
    dev.off()
}
cat(sprintf("[%s] written: %s  (total %.1f s)\n", TxID, output_path_plot,
            as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
