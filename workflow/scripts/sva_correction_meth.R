#!/usr/bin/env Rscript
# =============================================================================
# batch_correct_meth_fsva.R — LEAKAGE-SAFE methylation batch correction
#
# Frozen SVA (Parker, Corrada Bravo & Leek 2014, PeerJ 2:e561).
#
# STAGE 1 (subset, in memory):
#   sva()  on TRAIN samples, top-variable CpGs  -> train SVs
#   fsva() -> projects TEST samples onto the frozen train SV basis -> test SVs
#   Surrogate variables are SAMPLE-LEVEL vectors (length n_train / n_test),
#   so a subset of CpGs is sufficient to estimate them. A stability check
#   (--check_stability) verifies this empirically.
#
# STAGE 2 (genome-wide, chunked):
#   For EVERY CpG (all 18.7M), regress TRAIN M-values on [mod | train SVs],
#   take the SV coefficients (gammahat, per-CpG), and subtract:
#       train_corr = train - gammahat %*% t(train_SVs)
#       test_corr  = test  - gammahat %*% t(test_SVs)     <- SAME gammahat
#   No CpG is excluded from correction; the subset only fixed the SVs.
#
# NOTE ON CONFOUNDING: dataset (CPGEA/MCRPC) is aliased with metastasis.
# With mod = ~1 the SVs absorb the dataset axis, so removing them removes
# the metastasis biology too. N<->T (both CPGEA) is the clean contrast.
# =============================================================================
suppressPackageStartupMessages({
  library(sva); library(HDF5Array); library(rhdf5)
  library(DelayedMatrixStats); library(matrixStats)
  library(GenomicRanges); library(ggplot2)
})

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

in_h5        <- args$in_h5
in_rowranges <- args$in_rowranges
in_coldata   <- args$in_coldata
out_h5       <- args$out_h5
out_fit      <- args$out_fit
out_plot_ds  <- args$out_plot_ds
out_plot_cnd <- args$out_plot_cond
n_sv         <- if (!is.null(args$n_sv)) as.integer(args$n_sv) else 2L
n_top        <- if (!is.null(args$n_top)) as.integer(args$n_top) else 50000L
check_stab   <- !is.null(args$check_stability)

logf <- function(fmt, ...) {
  g <- gc(verbose = FALSE)
  cat(sprintf("[%s | %5.1f GB] %s\n", format(Sys.time(), "%H:%M:%S"),
              sum(g[, 2]) / 1024, sprintf(fmt, ...))); flush.console()
}

# ---- metadata ---------------------------------------------------------------
cd <- readRDS(in_coldata)
rr <- readRDS(in_rowranges)

samples   <- rownames(cd)
condition <- factor(cd$condition, levels = c("Normal", "Tumour", "metastasis"))
dataset   <- factor(cd$dataset,   levels = c("CPGEA", "MCRPC"))
split     <- factor(cd$split,     levels = c("train", "test"))
is_train  <- split == "train"

M_all <- HDF5Array::HDF5Array(in_h5, "M")
colnames(M_all) <- samples
stopifnot(ncol(M_all) == length(samples), nrow(M_all) == length(rr))

logf("M: %d CpGs x %d samples (train=%d, test=%d)",
     nrow(M_all), ncol(M_all), sum(is_train), sum(!is_train))

# =============================================================================
# STAGE 1 — estimate the surrogate variables (subset; SVs are sample-level)
# =============================================================================
logf("Computing per-CpG variance (lazy) to pick top %d ...", n_top)
rv <- DelayedMatrixStats::rowVars(M_all)

pick_top <- function(k) head(order(rv, decreasing = TRUE), k)

fit_sv <- function(idx, label) {
  Msub <- as.matrix(M_all[idx, ])
  mod  <- model.matrix(~ 1, data = data.frame(condition))   # protect nothing
  mod0 <- model.matrix(~ 1, data = data.frame(condition))
  sv   <- sva(Msub[, is_train, drop = FALSE], mod[is_train, , drop = FALSE],
              mod0[is_train, , drop = FALSE], n.sv = n_sv)
  logf("  [%s] SVA returned %d SV(s) from %d CpGs", label, sv$n.sv, length(idx))
  list(sv = sv, Msub = Msub, mod = mod)
}

# ---- optional: does the subset choice matter? -------------------------------
if (check_stab) {
  logf("Stability check: comparing SVs across subsets ...")
  set.seed(48047510)
  s_small <- fit_sv(pick_top(20000), "top20k")$sv$sv
  s_big   <- fit_sv(pick_top(200000), "top200k")$sv$sv
  s_rand  <- fit_sv(sort(sample(nrow(M_all), 50000)), "random50k")$sv$sv
  for (j in seq_len(n_sv)) {
    logf("  SV%d  cor(top20k, top200k) = %.4f | cor(random50k, top200k) = %.4f",
         j, abs(cor(s_small[, j], s_big[, j])), abs(cor(s_rand[, j], s_big[, j])))
  }
  rm(s_small, s_big, s_rand); invisible(gc(verbose = FALSE))
}

logf("Estimating SVs on TRAIN (top %d CpGs) ...", n_top)
top_idx <- pick_top(n_top)
res     <- fit_sv(top_idx, sprintf("top%dk", n_top / 1000))
svobj   <- res$sv
mod     <- res$mod

# ---- frozen SVA: project the TEST samples onto the train SV basis ----
logf("fsva (exact): projecting TEST samples onto the frozen train SVs ...")
fsvaobj <- fsva(dbdat  = res$Msub[, is_train, drop = FALSE],
                mod    = mod[is_train, , drop = FALSE],
                sv     = svobj,
                newdat = res$Msub[, !is_train, drop = FALSE],
                method = "exact")

SV_train <- svobj$sv          # n_train x n_sv
SV_test  <- fsvaobj$newsv     # n_test  x n_sv
stopifnot(nrow(SV_train) == sum(is_train), nrow(SV_test) == sum(!is_train))
logf("SVs: train %d x %d | test %d x %d",
     nrow(SV_train), ncol(SV_train), nrow(SV_test), ncol(SV_test))

rm(res); invisible(gc(verbose = FALSE))

# =============================================================================
# STAGE 2 — apply to ALL CpGs, chromosome by chromosome
# =============================================================================
chrom_per_cpg <- as.character(seqnames(rr))
chroms        <- unique(chrom_per_cpg)

if (file.exists(out_h5)) file.remove(out_h5)
rhdf5::h5createFile(out_h5)
rhdf5::h5createDataset(out_h5, "M", dims = c(nrow(M_all), ncol(M_all)),
                       chunk = c(min(10000, nrow(M_all)), ncol(M_all)),
                       level = 0, storage.mode = "double")

# design used to derive gammahat: [mod | SVs] on TRAIN, exactly as fsva does
D_train <- cbind(mod[is_train, , drop = FALSE], SV_train)
n_mod   <- ncol(mod)
sv_cols <- (n_mod + 1):(n_mod + ncol(SV_train))
XtX_inv <- solve(t(D_train) %*% D_train)

train_cols <- which(is_train)
test_cols  <- which(!is_train)

for (chrom in chroms) {
  idx <- which(chrom_per_cpg == chrom)
  stopifnot(identical(idx, seq(min(idx), max(idx))))   # contiguous in the H5

  logf("[%s] %d CpGs - correcting ...", chrom, length(idx))
  M <- as.matrix(M_all[idx, ])                          # CpGs x samples

  # per-CpG coefficients on the SVs, estimated from TRAIN only
  gammahat <- (M[, train_cols, drop = FALSE] %*% D_train %*% XtX_inv)[, sv_cols, drop = FALSE]

  # subtract the SV contribution: train with train SVs, test with test SVs
  M[, train_cols] <- M[, train_cols] - gammahat %*% t(SV_train)
  M[, test_cols]  <- M[, test_cols]  - gammahat %*% t(SV_test)

  rhdf5::h5write(M, out_h5, "M", index = list(idx, NULL))
  rm(M, gammahat); invisible(gc(verbose = FALSE))
  logf("[%s] done.", chrom)
}
rhdf5::H5close()

saveRDS(list(sv_train = SV_train, sv_test = SV_test, svobj = svobj,
             mod = mod, n_sv = n_sv, n_top = n_top, top_idx = top_idx,
             train_samples = samples[is_train], test_samples = samples[!is_train]),
        out_fit)

# =============================================================================
# QC PCA — fit on TRAIN, project TEST (before vs after)
# =============================================================================
logf("QC PCA ...")
M_corr <- HDF5Array::HDF5Array(out_h5, "M"); colnames(M_corr) <- samples

pca_pair <- function(Mmat, stage) {
  vg  <- head(order(matrixStats::rowVars(Mmat), decreasing = TRUE), 2000)
  Xtr <- t(Mmat[vg, is_train, drop = FALSE])
  ctr <- colMeans(Xtr)
  p   <- prcomp(Xtr, center = TRUE, scale. = FALSE)
  ve  <- round(100 * p$sdev^2 / sum(p$sdev^2), 1)
  Xte <- scale(t(Mmat[vg, !is_train, drop = FALSE]), center = ctr, scale = FALSE) %*% p$rotation
  rbind(
    data.frame(PC1 = p$x[, 1], PC2 = p$x[, 2], condition = condition[is_train],
               dataset = dataset[is_train], set = "train", stage = stage),
    data.frame(PC1 = Xte[, 1], PC2 = Xte[, 2], condition = condition[!is_train],
               dataset = dataset[!is_train], set = "test", stage = stage))
}

Mb <- as.matrix(M_all[top_idx, ]); Ma <- as.matrix(M_corr[top_idx, ])
df <- rbind(pca_pair(Mb, "before"), pca_pair(Ma, "after"))
df$stage <- factor(df$stage, levels = c("before", "after"))
rm(Mb, Ma); invisible(gc(verbose = FALSE))

mk <- function(col, vals, title) {
  ggplot(df, aes(PC1, PC2, color = .data[[col]])) +
    geom_point(size = 1.6, alpha = .8) +
    facet_grid(set ~ stage, scales = "free") +
    scale_color_manual(values = vals) + theme_bw(base_size = 12) +
    labs(title = title, subtitle = "PCA fit on train; test projected")
}
ggsave(out_plot_ds,
       mk("dataset", c(CPGEA = "#4575b4", MCRPC = "#d73027"),
          "Methylation fsva - by DATASET"), width = 10, height = 7)
ggsave(out_plot_cnd,
       mk("condition", c(Normal = "#2166ac", Tumour = "#d6604d", metastasis = "#1a9641"),
          "Methylation fsva - by CONDITION"), width = 10, height = 7)

logf("Saved: %s", out_h5)