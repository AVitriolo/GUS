# =============================================================================
# WGBS integration: CPGEA (Normal/Tumour) + MCRPC (metastasis)
#   LEAKAGE-SAFE: restricted to PAIRED samples (RNA & WGBS), with the
#   coverage filter and the imputation medians FIT ON TRAIN ONLY and
#   applied to test.
#
#   Per chromosome:
#     subset -> align/merge -> keep PAIRED samples ->
#     coverage filter (TRAIN) -> drop rows with no TRAIN median ->
#     impute with TRAIN per-condition medians -> save
# =============================================================================
args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

cpgea             <- args$cpgea
mcrpc             <- args$mcrpc
train_samples_path<- args$train_samples
test_samples_path <- args$test_samples
out_h5            <- args$out_h5
out_rowranges     <- args$out_rowranges
out_coldata       <- args$out_coldata
out_fit           <- args$out_fit          # stored train-derived fit

suppressPackageStartupMessages({
  library(SummarizedExperiment); library(GenomicRanges); library(HDF5Array)
  library(DelayedMatrixStats);   library(matrixStats);   library(dplyr)
  library(rhdf5)
})

if (requireNamespace("RhpcBLASctl", quietly = TRUE))
  RhpcBLASctl::blas_set_num_threads(8)

# ---- Config -----------------------------------------------------------------

PATHS <- list(local = list(cpgea = cpgea, mcrpc = mcrpc))
OUT_H5 <- out_h5; OUT_ROWRANGES <- out_rowranges
OUT_COLDATA <- out_coldata; OUT_FIT <- out_fit
CHECKPOINT_DIR <- file.path(dirname(out_h5), "checkpoints_median_percond")

FILTER_CFG <- list(
  Normal     = list(min_cov = 10, min_frac = 0.80),
  Tumour     = list(min_cov = 10, min_frac = 0.80),
  metastasis = list(min_cov = 5,  min_frac = 0.80)
)

TRAIN_SAMPLES <- readLines(train_samples_path)
TEST_SAMPLES  <- readLines(test_samples_path)
PAIRED        <- c(TRAIN_SAMPLES, TEST_SAMPLES)

# ---- Logging ----------------------------------------------------------------

.PIPE_T0 <- NULL
log_start <- function() {
  .PIPE_T0 <<- Sys.time()
  logf("Pipeline started at %s", format(.PIPE_T0, "%Y-%m-%d %H:%M:%S"))
}
logf <- function(fmt, ...) {
  msg <- sprintf(fmt, ...)
  g <- gc(verbose = FALSE)
  cat(sprintf("[%s | %5.1f GB] %s\n", format(Sys.time(), "%H:%M:%S"),
              sum(g[, 2]) / 1024, msg)); flush.console()
}
log_step <- function(name) { cat("\n"); logf("==== %s ====", name) }

check_mask <- function(mask, name = "mask") {
  na <- sum(is.na(mask))
  logf("  mask '%s': length=%d, TRUE=%d, FALSE=%d, NA=%d",
       name, length(mask), sum(mask, na.rm = TRUE), sum(!mask, na.rm = TRUE), na)
  if (na > 0) stop(sprintf("Mask '%s' contains %d NAs.", name, na))
  invisible(TRUE)
}

ensure_dirs <- function() {
  for (d in c(CHECKPOINT_DIR, dirname(OUT_H5), dirname(OUT_FIT)))
    if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ---- Helpers ----------------------------------------------------------------

coord_key <- function(se) {
  rr <- rowRanges(se)
  paste0(as.character(seqnames(rr)), "_", start(rr), "_", end(rr))
}

load_datasets <- function(paths) list(
  CPGEA = HDF5Array::loadHDF5SummarizedExperiment(paths$cpgea),
  MCRPC = HDF5Array::loadHDF5SummarizedExperiment(paths$mcrpc)
)

subset_to_chromosome <- function(se, chrom)
  se[which(as.character(seqnames(rowRanges(se))) == chrom), ]

align_and_merge <- function(se_a, se_b, name_a = "CPGEA", name_b = "MCRPC") {
  key_a <- coord_key(se_a); key_b <- coord_key(se_b)
  stopifnot(!anyDuplicated(key_a), !anyDuplicated(key_b))
  common <- intersect(key_a, key_b)
  se_a2 <- se_a[match(common, key_a), ]; se_b2 <- se_b[match(common, key_b), ]
  rownames(se_a2) <- common; rownames(se_b2) <- common
  merged <- combineCols(se_a2, se_b2, use.names = TRUE)
  colData(merged)$dataset <- factor(
    c(rep(name_a, ncol(se_a2)), rep(name_b, ncol(se_b2))),
    levels = c(name_a, name_b))
  merged
}

infer_condition <- function(se) {
  nms <- rownames(colData(se))
  colData(se)$condition <- factor(dplyr::case_when(
    grepl("^N", nms)   ~ "Normal",
    grepl("^T", nms)   ~ "Tumour",
    grepl("^DTB", nms) ~ "metastasis",
    TRUE               ~ NA_character_
  ), levels = c("Normal", "Tumour", "metastasis"))
  se
}

#' Keep only the PAIRED samples (RNA & WGBS), and tag train/test.
subset_to_paired <- function(se) {
  keep <- colnames(se) %in% PAIRED
  se <- se[, keep]
  colData(se)$split <- factor(
    ifelse(colnames(se) %in% TRAIN_SAMPLES, "train", "test"),
    levels = c("train", "test"))
  se
}

# ---- [TRAIN-ONLY] coverage filter -------------------------------------------
#' The CpG universe is decided using TRAIN samples only; the resulting
#' mask is then applied to every sample (train and test alike).
rows_passing_coverage_train <- function(se, cov_assay = "cov", cfg = FILTER_CFG) {
  cov   <- assay(se, cov_assay)
  cond  <- colData(se)$condition
  split <- colData(se)$split
  pass  <- rep(TRUE, nrow(se))

  for (cnd in names(cfg)) {
    cols <- which(cond == cnd & split == "train")     # <-- TRAIN ONLY
    if (!length(cols)) {
      warning(sprintf("No TRAIN samples for condition '%s'; skipping.", cnd)); next
    }
    frac_ok <- DelayedMatrixStats::rowMeans2(
      cov[, cols, drop = FALSE] >= cfg[[cnd]]$min_cov, na.rm = TRUE)
    this_pass <- !is.nan(frac_ok) & (frac_ok >= cfg[[cnd]]$min_frac)
    logf("    %-10s: %d / %d CpGs pass (>= %.0f%% of %d TRAIN samples, cov >= %dx)",
         cnd, sum(this_pass), length(this_pass),
         100 * cfg[[cnd]]$min_frac, length(cols), cfg[[cnd]]$min_cov)
    pass <- pass & this_pass
  }
  check_mask(pass, "coverage_pass_train")
  pass
}

# ---- [TRAIN-ONLY] rows that have a usable train median -----------------------
#' Imputation uses TRAIN per-condition medians. A CpG is only usable if
#' every condition has >= 1 non-NA TRAIN value; otherwise   no median exists.
rows_with_train_median <- function(se, assay_name = "beta") {
  mat   <- assay(se, assay_name)
  cond  <- colData(se)$condition
  split <- colData(se)$split
  keep  <- rep(TRUE, nrow(mat))

  for (cnd in levels(cond)) {
    cols <- which(cond == cnd & split == "train")
    if (!length(cols)) next
    n_obs <- DelayedMatrixStats::rowSums2(!is.na(mat[, cols, drop = FALSE]))
    no_median <- n_obs == 0
    logf("    no TRAIN median for %-10s: %d rows", cnd, sum(no_median))
    keep <- keep & !no_median
  }
  check_mask(keep, "has_train_median")
  keep
}

beta_to_m <- function(beta, eps = 1e-6) {
  beta <- pmin(pmax(beta, eps), 1 - eps)
  log2(beta / (1 - beta))
}

# ---- [TRAIN-FIT / APPLY-TO-ALL] imputation ----------------------------------
#' Medians are computed from TRAIN samples only, then used to fill NAs in
#' BOTH train and test. Test values never influence the fill values.
impute_with_train_medians <- function(se_chr, chrom, assay_name = "beta") {
  logf("[%s] Realising %d x %d beta (~%.2f GB)...", chrom,
       nrow(se_chr), ncol(se_chr),
       as.numeric(nrow(se_chr)) * ncol(se_chr) * 8 / 1e9)
  t0   <- Sys.time()
  beta <- as.matrix(assay(se_chr, assay_name))
  logf("[%s]   realised in %.1f s", chrom,
       as.numeric(difftime(Sys.time(), t0, units = "secs")))

  cond  <- as.character(colData(se_chr)$condition)
  split <- as.character(colData(se_chr)$split)

  logf("[%s] NAs to impute: %d (%.2f%%)", chrom, sum(is.na(beta)),
       100 * sum(is.na(beta)) / length(beta))

  medians <- list()
  t0 <- Sys.time()
  for (cnd in unique(cond)) {
    tr_cols  <- which(cond == cnd & split == "train")
    all_cols <- which(cond == cnd)

    # ---- FIT: median from TRAIN samples of this condition ----
    row_med <- matrixStats::rowMedians(beta[, tr_cols, drop = FALSE], na.rm = TRUE)
    medians[[cnd]] <- row_med

    # ---- APPLY: fill NAs in ALL samples of this condition (train + test) ----
    sub    <- beta[, all_cols, drop = FALSE]
    na_idx <- which(is.na(sub), arr.ind = TRUE)
    if (nrow(na_idx) > 0) {
      sub[na_idx]      <- row_med[na_idx[, "row"]]
      beta[, all_cols] <- sub
      logf("[%s]   %s: filled %d NAs using medians from %d TRAIN samples",
           chrom, cnd, nrow(na_idx), length(tr_cols))
    }
  }
  logf("[%s] Imputation done in %.1f s", chrom,
       as.numeric(difftime(Sys.time(), t0, units = "secs")))

  if (sum(is.na(beta)) > 0)
    stop(sprintf("[%s] %d NAs remain after train-median imputation.",
                 chrom, sum(is.na(beta))))

  list(beta = beta, medians = medians)
}

# ---- One chromosome ---------------------------------------------------------

process_chromosome <- function(chrom, paths, dir = CHECKPOINT_DIR, cfg = FILTER_CFG) {
  out_file <- file.path(dir, sprintf("chr_%s.rds", chrom))
  if (file.exists(out_file)) {
    logf("[%s] Already done -- skipping", chrom); return(invisible(NULL))
  }
  log_step(sprintf("Chromosome %s", chrom))

  ds <- load_datasets(paths)
  cpgea_chr <- subset_to_chromosome(ds$CPGEA, chrom)
  mcrpc_chr <- subset_to_chromosome(ds$MCRPC, chrom)
  logf("[%s] CPGEA=%d CpGs, MCRPC=%d CpGs", chrom, nrow(cpgea_chr), nrow(mcrpc_chr))
  rm(ds); invisible(gc(verbose = FALSE))

  merged <- align_and_merge(cpgea_chr, mcrpc_chr, "CPGEA", "MCRPC")
  logf("[%s] Common CpGs: %d", chrom, nrow(merged))
  rm(cpgea_chr, mcrpc_chr); invisible(gc(verbose = FALSE))

  merged <- infer_condition(merged)

  # ---- keep PAIRED samples only, tag train/test ----
  n_before <- ncol(merged)
  merged   <- subset_to_paired(merged)
  logf("[%s] Samples: %d -> %d paired (train=%d, test=%d)", chrom,
       n_before, ncol(merged),
       sum(colData(merged)$split == "train"),
       sum(colData(merged)$split == "test"))
  stopifnot(ncol(merged) == length(PAIRED))

  logf("[%s] Coverage filter (TRAIN only):", chrom)
  merged <- merged[rows_passing_coverage_train(merged, "cov", cfg), ]
  logf("[%s] Post-coverage: %d CpGs", chrom, nrow(merged))

  logf("[%s] Dropping rows without a TRAIN median:", chrom)
  merged <- merged[rows_with_train_median(merged, "beta"), ]
  logf("[%s] Post-median-check: %d CpGs", chrom, nrow(merged))

  logf("[%s] Imputation (TRAIN per-condition medians):", chrom)
  imp <- impute_with_train_medians(merged, chrom)

  saveRDS(list(chrom    = chrom,
               row_keys = rownames(merged),
               colData  = as.data.frame(colData(merged)),
               beta_imp = imp$beta,
               medians  = imp$medians),      # train-derived fit, auditable
          file = out_file)

  rm(merged, imp); invisible(gc(verbose = FALSE))
  logf("[%s] Done.", chrom)
}

# ---- Reassembly -------------------------------------------------------------

keys_to_granges <- function(keys) {
  m <- do.call(rbind, strsplit(keys, "_"))
  GRanges(seqnames = m[, 1],
          ranges = IRanges(start = as.integer(m[, 2]), end = as.integer(m[, 3])))
}

reassemble_to_hdf5 <- function(dir = CHECKPOINT_DIR, out_h5 = OUT_H5,
                               out_rr = OUT_ROWRANGES, out_cd = OUT_COLDATA,
                               out_fit = OUT_FIT,
                               chr_order = paste0("chr", c(1:22, "X", "Y"))) {
  files <- file.path(dir, sprintf("chr_%s.rds", chr_order))
  if (any(!file.exists(files)))
    stop("Missing checkpoints: ", paste(chr_order[!file.exists(files)], collapse = ", "))

  obj1     <- readRDS(files[1])
  n_cols   <- nrow(obj1$colData)
  col_data <- obj1$colData
  rm(obj1); invisible(gc(verbose = FALSE))

  n_rows <- sum(sapply(files, function(f) length(readRDS(f)$row_keys)))
  logf("Total: %d CpGs x %d samples (~%.1f GB per assay)", n_rows, n_cols,
       as.numeric(n_rows) * n_cols * 8 / 1e9)

  if (file.exists(out_h5)) file.remove(out_h5)
  rhdf5::h5createFile(out_h5)
  for (nm in c("beta", "M"))
    rhdf5::h5createDataset(out_h5, nm, dims = c(n_rows, n_cols),
                           chunk = c(1000, n_cols), level = 0,
                           storage.mode = "double")

  row_offset <- 0L; all_keys <- character(0); fit_medians <- list()

  for (f in files) {
    chrom <- sub("chr_(chr.+)\\.rds", "\\1", basename(f))
    obj <- readRDS(f)
    nr  <- nrow(obj$beta_imp)
    idx <- seq(row_offset + 1L, row_offset + nr)

    rhdf5::h5write(obj$beta_imp,            out_h5, "beta", index = list(idx, NULL))
    rhdf5::h5write(beta_to_m(obj$beta_imp), out_h5, "M",    index = list(idx, NULL))

    all_keys           <- c(all_keys, obj$row_keys)
    fit_medians[[chrom]] <- obj$medians
    row_offset         <- row_offset + nr
    rm(obj); invisible(gc(verbose = FALSE))
    logf("[%s] written (offset %d)", chrom, row_offset)
  }
  rhdf5::H5close()

  saveRDS(keys_to_granges(all_keys), out_rr)
  saveRDS(col_data, out_cd)
  saveRDS(list(train_samples = TRAIN_SAMPLES,
               test_samples  = TEST_SAMPLES,
               filter_cfg    = FILTER_CFG,
               medians       = fit_medians),
          out_fit)
  logf("Metadata + train-fit saved.")
}

load_integrated_se <- function(out_h5 = OUT_H5, out_rr = OUT_ROWRANGES,
                               out_cd = OUT_COLDATA) {
  SummarizedExperiment(
    assays    = list(beta = HDF5Array::HDF5Array(out_h5, "beta"),
                     M    = HDF5Array::HDF5Array(out_h5, "M")),
    rowRanges = readRDS(out_rr),
    colData   = S4Vectors::DataFrame(readRDS(out_cd)))
}

# ---- Driver -----------------------------------------------------------------

run_pipeline <- function(paths, only_chrom = NULL) {
  log_start(); ensure_dirs()
  logf("Paired samples: %d (train=%d, test=%d)",
       length(PAIRED), length(TRAIN_SAMPLES), length(TEST_SAMPLES))

  if (identical(only_chrom, "REASSEMBLE")) {
    log_step("Reassemble to HDF5"); reassemble_to_hdf5()
    logf("Done in %.1f min",
         as.numeric(difftime(Sys.time(), .PIPE_T0, units = "mins")))
    return(invisible(NULL))
  }

  chrom_list <- if (!is.null(only_chrom)) only_chrom
                else paste0("chr", c(2, 1, 3:22, "X", "Y"))
  for (chrom in chrom_list)
    process_chromosome(chrom, paths, CHECKPOINT_DIR, FILTER_CFG)

  if (is.null(only_chrom)) {
    log_step("Reassemble to HDF5"); reassemble_to_hdf5()
  }
  logf("Done in %.1f min",
       as.numeric(difftime(Sys.time(), .PIPE_T0, units = "mins")))
  invisible(NULL)
}

run_pipeline(PATHS$local)