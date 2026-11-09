#!/usr/bin/env Rscript
# =============================================================================
# build_purebeta_md_store.R
# Three MethylDriver input stores on NON-ComBat betas, same CpGs, same samples, same missing cells:
#   bulk     : v2 raw beta (integrated_wgbs.h5)
#   camdac   : normals v2 raw; tumours/metastases CAMDAC m_t (Prostate_GAS_purified integrated_wgbs.h5, pre-ComBat beta)
#   purebeta : normals v2 raw; tumours/metastases PureBeta on v2 raw, tumours and metastases corrected as two separate
#              cohorts (PAMES purity, NT / M panels). PureBeta corrects each CpG independently; values are not clipped
# CpGs : the BED's CpGs present in both v2 and the CAMDAC (v1) store. Samples: present in both stores.
# Mask : cells imputed in either store (value = the condition's train median, how both builds imputed) -> NA in all three;
#        PureBeta is run on the masked raw values.
# Writes <out_dir>/{bulk,camdac,purebeta}.h5 (dataset "M" = beta, CpGs x samples) + _rowRanges.rds + _colData.rds
# Usage: Rscript build_purebeta_md_store.R --bed=<percpg.bed> --out_dir=<dir> [--cores=100]
# =============================================================================
suppressPackageStartupMessages({ library(PureBeta); library(HDF5Array); library(GenomicRanges); library(data.table) })
args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
cores <- as.integer(if (is.null(args$cores)) 100 else args$cores); out <- args$out_dir; dir.create(out, recursive = TRUE, showWarnings = FALSE)
t0 <- Sys.time(); tm <- function() sprintf("%.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))
V <- "resources/integrated_wgbs_merged_v2_hdf5/integrated_wgbs"; C <- "/home/ieo7293/Prostate_GAS_purified/resources/integrated_wgbs_hdf5/integrated_wgbs"
rrV <- readRDS(paste0(V, "_rowRanges.rds")); cdV <- readRDS(paste0(V, "_colData.rds")); idV <- paste0(seqnames(rrV), "_", start(rrV))
rrC <- readRDS(paste0(C, "_rowRanges.rds")); cdC <- readRDS(paste0(C, "_colData.rds")); idC <- paste0(seqnames(rrC), "_", start(rrC))
want <- unique(fread(args$bed, header = FALSE)$V4); want <- want[want %in% idV & want %in% idC]
iV <- match(want, idV); o <- order(iV); want <- want[o]; iV <- iV[o]; iC <- match(want, idC)
smp <- intersect(rownames(cdV), rownames(cdC)); cd <- cdV[smp, ]
cat(sprintf("CpGs %d | samples %d (%s)\n", length(want), length(smp), paste(names(table(cd$condition)), table(cd$condition), collapse = " ")))
rd <- function(f, ds, i, s, cols) { X <- as.matrix(HDF5Array(paste0(f, ".h5"), ds)[i, match(s, cols)]); dimnames(X) <- list(want, s); X }
RAW <- rd(V, "beta", iV, smp, rownames(cdV)); CAM <- rd(C, "beta", iC, smp, rownames(cdC)); cat("read,", tm(), "\n")

# ---- imputed cells: equal to the condition's train median ----
imputed <- function(X, cdx) {
    m <- matrix(FALSE, nrow(X), ncol(X))
    for (k in unique(cdx$condition)) {
        j <- which(cdx[colnames(X), "condition"] == k); tr <- j[cdx[colnames(X)[j], "split"] == "train"]
        med <- apply(X[, tr, drop = FALSE], 1, median, na.rm = TRUE); m[, j] <- abs(X[, j] - med) < 1e-6
    }
    m[is.na(m)] <- FALSE; m
}
mV <- imputed(RAW, as.data.frame(cdV)); mC <- imputed(CAM, as.data.frame(cdC))
cond <- as.character(cd$condition); nn <- cond == "Normal"
mC[, nn] <- FALSE                                     # CAMDAC store normals are not used
mask <- mV | mC | is.na(RAW) | (is.na(CAM) & matrix(!nn, nrow(CAM), ncol(CAM), byrow = TRUE))
cat(sprintf("masked cells: imputed v2 %.2f%% | imputed CAMDAC (T/M) %.2f%% | total %.2f%%\n", 100 * mean(mV), 100 * mean(mC[, !nn]), 100 * mean(mask)))
RAW[mask] <- NA
CAMs <- RAW; CAMs[, !nn] <- CAM[, !nn]; CAMs[mask] <- NA

# ---- PureBeta, tumours and metastases as separate cohorts ----
pn <- read.delim("resources/purity/pames_purity_NT.tsv"); pm <- read.delim("resources/purity/pames_purity_M.tsv")
pur <- c(setNames(pn$Purity, pn$Sample), setNames(pm[[2]], pm[[1]]))
# PureBeta's own per-CpG function (adjustBeta), called per CpG as beta_correction_for_cohorts does (same seed). A CpG
# whose flexmix mixture does not converge gets PureBeta's one-group formula directly (adjustBeta itself cannot run with
# one component: getModel() fails on a plain flexmix fit): lm(beta ~ 1 - purity), tumour = intercept + residual.
# Fewer than 5 values -> NA. Work is handed out in small chunks so slow CpGs do not hold up a whole prescheduled block.
one_group <- function(v, p, s) {
    y <- setNames(rep(NA_real_, length(s)), s); ok <- !is.na(v) & !is.na(p)
    if (sum(ok) < 5 || var(v[ok]) == 0) return(y)
    f <- lm(v[ok] ~ I(1 - p[ok])); y[ok] <- coef(f)[1] + residuals(f); y
}
adjust_one <- function(v, p, s, seed) {
    na <- list(y = setNames(rep(NA_real_, length(s)), s), st = "na")
    if (sum(!is.na(v)) < 5) return(na)
    r <- tryCatch(PureBeta:::adjustBeta(c(seed, v), p, s, seed = TRUE)$y.tum, error = function(e) NULL)
    if (!is.null(r)) return(list(y = r[s], st = "mixture"))
    tryCatch(list(y = one_group(v, p, s), st = "one group"), error = function(e) na)
}
# PSOCK workers (not forks): a forked R worker ends up copying the parent's ~3 GB of matrices (the garbage collector
# touches them), so 100 forks exhausted the machine's memory. Each PSOCK worker only receives its own chunk.
cl <- parallel::makePSOCKcluster(cores)
invisible(parallel::clusterEvalQ(cl, suppressPackageStartupMessages({ library(flexmix); library(PureBeta) })))
parallel::clusterExport(cl, c("adjust_one", "one_group"))
PB <- RAW; PB[, !nn] <- NA
for (k in c("Tumour", "metastasis")) {
    s <- smp[cond == k]; s <- s[s %in% names(pur) & !is.na(pur[s])]; X <- RAW[, s]; seed <- 2000 + nrow(X); p <- pur[s]
    t1 <- Sys.time(); ch <- lapply(split(seq_len(nrow(X)), ceiling(seq_len(nrow(X)) / 200)), function(jj) X[jj, , drop = FALSE])
    res <- parallel::parLapplyLB(cl, ch, function(x, p, s, seed) lapply(seq_len(nrow(x)), function(j) suppressWarnings(adjust_one(x[j, ], p, s, seed))),
                                 p = p, s = s, seed = seed)
    res <- unlist(res, recursive = FALSE); PB[, s] <- do.call(rbind, lapply(res, `[[`, "y"))
    st <- table(vapply(res, `[[`, character(1), "st"))
    cat(sprintf("PureBeta %s: %d samples (%d without purity -> NA), %.1f min | %s\n", k, length(s), sum(cond == k) - length(s),
                as.numeric(difftime(Sys.time(), t1, units = "mins")), paste(names(st), st, collapse = ", ")))
}
PB[mask] <- NA

# ---- write ----
for (nm in c("bulk", "camdac", "purebeta")) {
    X <- switch(nm, bulk = RAW, camdac = CAMs, purebeta = PB); f <- file.path(out, paste0(nm, ".h5")); if (file.exists(f)) file.remove(f)
    writeHDF5Array(X, filepath = f, name = "M", chunkdim = c(min(1000, nrow(X)), ncol(X)))
    saveRDS(rrV[iV], sub("\\.h5$", "_rowRanges.rds", f)); saveRDS(cd, sub("\\.h5$", "_colData.rds", f))
    cat(sprintf("%-9s NA %.2f%% | T/M values < 0 or > 1: %.2f%% | at exactly 0 or 1: %.2f%%\n", nm, 100 * mean(is.na(X)),
                100 * mean(X[, !nn] < 0 | X[, !nn] > 1, na.rm = TRUE), 100 * mean(X[, !nn] == 0 | X[, !nn] == 1, na.rm = TRUE)))
}
parallel::stopCluster(cl); cat("done,", tm(), "\n")
