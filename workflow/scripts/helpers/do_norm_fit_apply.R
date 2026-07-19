donorm_fit <- function(train_counts, method = "TMM", ref_name = NULL) {

    # choose the reference explicitly and deterministically.
    # default: the train sample whose 75th percentile (of non-zero counts)
    # is closest to the mean 75th percentile across train samples.
    if (is.null(ref_name)) {
        f75 <- apply(train_counts, 2, function(x) quantile(x[x > 0], 0.75))
        ref_idx  <- which.min(abs(f75 - mean(f75)))
        ref_name <- colnames(train_counts)[ref_idx]
    } else {
        stopifnot(ref_name %in% colnames(train_counts))
        ref_idx <- match(ref_name, colnames(train_counts))
    }

    d <- edgeR::DGEList(counts = train_counts)
    d <- edgeR::calcNormFactors(d, method = method, refColumn = ref_idx)  # EXPLICIT

    e <- edgeR::cpm(d, normalized.lib.sizes = TRUE)
    libsize_const <- mean(d$samples$lib.size)
    e <- e * libsize_const / 1e6

    list(
        normalized    = as.data.frame(e),
        ref_profile   = train_counts[, ref_idx, drop = FALSE],  # the raw reference column
        ref_name      = ref_name,
        libsize_const = libsize_const,
        method        = method
    )
}

donorm_apply <- function(test_counts, fit) {

    # prepend the train reference profile so refColumn = 1 is unambiguous
    combined <- cbind(fit$ref_profile, test_counts)
    colnames(combined)[1] <- paste0("__REF__", fit$ref_name)

    d <- edgeR::DGEList(counts = combined)
    d <- edgeR::calcNormFactors(d, method = fit$method, refColumn = 1)   # EXPLICIT: train ref

    e <- edgeR::cpm(d, normalized.lib.sizes = TRUE)
    e <- e * fit$libsize_const / 1e6      # SAME constant as train, not test's own

    e <- e[, -1, drop = FALSE]            # drop the reference column
    as.data.frame(e)
}
