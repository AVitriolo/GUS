# =============================================================================
# WGBS integration: CPGEA (Normal/Tumour) + MCRPC (metastasis)
#   Per-chromosome pipeline. Each chromosome is fully self-contained:
#     subset originals -> align -> per-condition coverage filter ->
#     drop all-NA rows -> methyLImp2 impute -> save.
#   No global filter step. Each chromosome's output is self-describing
#   (row keys + colData + imputed beta), so chromosomes are independent
#   and the pipeline is naturally parallel for SLURM arrays.
# =============================================================================
args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                     #  read args

cpgea               <- args$cpgea
mcrpc               <- args$mcrpc
out_h5              <- args$out_h5
out_rowranges       <- args$out_rowranges
out_coldata         <- args$out_coldata

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(GenomicRanges)
  library(HDF5Array)
  library(DelayedMatrixStats)
  library(matrixStats)
  library(dplyr)
  library(rhdf5)
})

options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!requireNamespace("ggVennDiagram", quietly = TRUE)) {
    install.packages("ggVennDiagram", repos = "https://cloud.r-project.org")
}

if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(8)
}

# ---- Config -----------------------------------------------------------------

PATHS <- list(local = list(cpgea = cpgea, mcrpc = mcrpc))
OUT_H5         <- out_h5
OUT_ROWRANGES  <- out_rowranges
OUT_COLDATA    <- out_coldata
CHECKPOINT_DIR <- file.path(dirname(out_h5), "checkpoints_median_percond")

FILTER_CFG <- list(
  Normal     = list(min_cov = 10, min_frac = 0.80),
  Tumour     = list(min_cov = 10, min_frac = 0.80),
  metastasis = list(min_cov = 5,  min_frac = 0.80)
)

# ---- Logging ----------------------------------------------------------------

.PIPE_T0 <- NULL

log_start <- function() {
  .PIPE_T0 <<- Sys.time()
  logf("Pipeline started at %s", format(.PIPE_T0, "%Y-%m-%d %H:%M:%S"))
}

logf <- function(fmt, ...) {
  msg <- sprintf(fmt, ...)
  g <- gc(verbose = FALSE)
  mem <- sprintf("%5.1f GB", sum(g[, 2]) / 1024)
  cat(sprintf("[%s | %s] %s\n", format(Sys.time(), "%H:%M:%S"), mem, msg))
  flush.console()
}

log_step <- function(name) {
  cat("\n"); logf("==== %s ====", name)
}

check_mask <- function(mask, name = "mask") {
  na <- sum(is.na(mask))
  logf("  mask '%s': length=%d, TRUE=%d, FALSE=%d, NA=%d",
       name, length(mask), sum(mask, na.rm = TRUE), sum(!mask, na.rm = TRUE), na)
  if (na > 0) stop(sprintf("Mask '%s' contains %d NAs.", name, na))
  invisible(TRUE)
}

ensure_dirs <- function() {
  if (!dir.exists(CHECKPOINT_DIR))
    dir.create(CHECKPOINT_DIR, recursive = TRUE, showWarnings = FALSE)
  out_dir <- dirname(OUT_H5)
  if (!dir.exists(out_dir))
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

# ---- Helpers ----------------------------------------------------------------

coord_key <- function(se) {
  rr <- rowRanges(se)
  paste0(as.character(seqnames(rr)), "_", start(rr), "_", end(rr))
}

load_datasets <- function(paths) {
  cpgea <- HDF5Array::loadHDF5SummarizedExperiment(paths$cpgea)
  mcrpc <- HDF5Array::loadHDF5SummarizedExperiment(paths$mcrpc)
  list(CPGEA = cpgea, MCRPC = mcrpc)
}

subset_to_chromosome <- function(se, chrom) {
  idx <- which(as.character(seqnames(rowRanges(se))) == chrom)
  se[idx, ]
}

align_and_merge <- function(se_a, se_b, name_a = "CPGEA", name_b = "MCRPC") {
  key_a <- coord_key(se_a); key_b <- coord_key(se_b)
  stopifnot(!anyDuplicated(key_a), !anyDuplicated(key_b))
  common <- intersect(key_a, key_b)
  se_a2 <- se_a[match(common, key_a), ]
  se_b2 <- se_b[match(common, key_b), ]
  rownames(se_a2) <- common; rownames(se_b2) <- common
  merged <- combineCols(se_a2, se_b2, use.names = TRUE)
  colData(merged)$dataset <- factor(
    c(rep(name_a, ncol(se_a2)), rep(name_b, ncol(se_b2))),
    levels = c(name_a, name_b)
  )
  merged
}

infer_condition <- function(se) {
  cd  <- as.data.frame(colData(se))
  nms <- rownames(cd)
  guess <- dplyr::case_when(
    grepl("^N", nms) ~ "Normal",
    grepl("^T", nms) ~ "Tumour",
    grepl("^D", nms) ~ "metastasis",
    TRUE             ~ NA_character_
  )
  current <- if ("condition" %in% colnames(cd)) as.character(cd$condition) else rep(NA_character_, nrow(cd))
  current[is.na(current)] <- guess[is.na(current)]
  colData(se)$condition <- factor(current, levels = c("Normal", "Tumour", "metastasis"))
  se
}

rows_passing_coverage <- function(se, cov_assay = "cov", cfg = FILTER_CFG,
                                  condition_col = "condition") {
  cov  <- assay(se, cov_assay)
  cond <- colData(se)[[condition_col]]
  pass <- rep(TRUE, nrow(se))
  for (cnd in names(cfg)) {
    cols <- which(cond == cnd)
    if (!length(cols)) {
      warning(sprintf("No samples for condition '%s'; skipping.", cnd)); next
    }
    frac_ok <- DelayedMatrixStats::rowMeans2(
      cov[, cols, drop = FALSE] >= cfg[[cnd]]$min_cov, na.rm = TRUE
    )
    this_pass <- !is.nan(frac_ok) & (frac_ok >= cfg[[cnd]]$min_frac)
    logf("    %-10s: %d / %d CpGs pass (>= %.0f%% samples with cov >= %dx)",
         cnd, sum(this_pass), length(this_pass),
         100 * cfg[[cnd]]$min_frac, cfg[[cnd]]$min_cov)
    pass <- pass & this_pass
  }
  check_mask(pass, "coverage_pass")
  pass
}

drop_all_na_in_any_dataset <- function(se, dataset_col = "dataset",
                                       assay_name = "beta") {
  mat <- assay(se, assay_name)
  ds  <- colData(se)[[dataset_col]]
  keep <- rep(TRUE, nrow(mat))
  for (d in unique(ds)) {
    cols <- which(ds == d)
    n_obs <- DelayedMatrixStats::rowSums2(!is.na(mat[, cols, drop = FALSE]))
    all_na <- n_obs == 0
    logf("    rows all-NA in %-7s: %d", d, sum(all_na))
    keep <- keep & !all_na
  }
  check_mask(keep, "drop_all_na_keep")
  se[keep, ]
}

beta_to_m <- function(beta, eps = 1e-6) {
  beta <- pmin(pmax(beta, eps), 1 - eps)
  log2(beta / (1 - beta))
}

impute_median_per_condition <- function(se_chr, chrom,
                                        condition_col = "condition",
                                        assay_name    = "beta") {
  logf("[%s] Realising %d x %d beta to dense (~%.2f GB)...",
       chrom, nrow(se_chr), ncol(se_chr),
       as.numeric(nrow(se_chr)) * ncol(se_chr) * 8 / 1e9)
  t0 <- Sys.time()
  beta <- as.matrix(assay(se_chr, assay_name))
  logf("[%s]   realised in %.1f s",
       chrom, as.numeric(difftime(Sys.time(), t0, units = "secs")))

  cond <- as.character(colData(se_chr)[[condition_col]])
  n_na_total <- sum(is.na(beta))
  logf("[%s] NAs to impute: %d (%.2f%%)",
       chrom, n_na_total,
       100 * n_na_total / (nrow(beta) * ncol(beta)))

  t0 <- Sys.time()
  for (cnd in unique(cond)) {
    cols <- which(cond == cnd)
    sub  <- beta[, cols, drop = FALSE]
    row_med <- matrixStats::rowMedians(sub, na.rm = TRUE)
    na_idx <- which(is.na(sub), arr.ind = TRUE)
    if (nrow(na_idx) > 0) {
      sub[na_idx] <- row_med[na_idx[, "row"]]
      beta[, cols] <- sub
      logf("[%s]   %s: filled %d NAs with per-CpG median (of %d samples)",
           chrom, cnd, nrow(na_idx), length(cols))
    }
  }
  logf("[%s] Imputation done in %.1f s",
       chrom, as.numeric(difftime(Sys.time(), t0, units = "secs")))

  remaining <- sum(is.na(beta))
  if (remaining > 0) {
    logf("[%s] WARNING: %d NAs remain.", chrom, remaining)
    row_med_global <- matrixStats::rowMedians(beta, na.rm = TRUE)
    na_idx <- which(is.na(beta), arr.ind = TRUE)
    beta[na_idx] <- row_med_global[na_idx[, "row"]]
    logf("[%s]   filled remaining with global per-row median; final NAs: %d",
         chrom, sum(is.na(beta)))
  }

  beta
}

# ---- One chromosome ---------------------------------------------------------

process_chromosome <- function(chrom, paths, dir = CHECKPOINT_DIR,
                               cfg = FILTER_CFG) {
  out_file <- file.path(dir, sprintf("chr_%s.rds", chrom))
  if (file.exists(out_file)) {
    logf("[%s] Already done at %s -- skipping", chrom, out_file)
    return(invisible(NULL))
  }

  log_step(sprintf("Chromosome %s", chrom))

  logf("[%s] Loading originals (lazy HDF5)...", chrom)
  ds <- load_datasets(paths)

  logf("[%s] Subsetting both SEs to %s...", chrom, chrom)
  cpgea_chr <- subset_to_chromosome(ds$CPGEA, chrom)
  mcrpc_chr <- subset_to_chromosome(ds$MCRPC, chrom)
  logf("[%s] CPGEA=%d CpGs, MCRPC=%d CpGs", chrom, nrow(cpgea_chr), nrow(mcrpc_chr))
  rm(ds); invisible(gc(verbose = FALSE))

  logf("[%s] Aligning & merging on common CpGs...", chrom)
  t0 <- Sys.time()
  merged <- align_and_merge(cpgea_chr, mcrpc_chr, "CPGEA", "MCRPC")
  logf("[%s]   common: %d CpGs (in %.1f s)",
       chrom, nrow(merged),
       as.numeric(difftime(Sys.time(), t0, units = "secs")))
  rm(cpgea_chr, mcrpc_chr); invisible(gc(verbose = FALSE))

  merged <- infer_condition(merged)

  logf("[%s] Per-condition coverage filter...", chrom)
  keep <- rows_passing_coverage(merged, cov_assay = "cov", cfg = cfg)
  merged <- merged[keep, ]
  logf("[%s] Post-coverage: %d CpGs", chrom, nrow(merged))

  logf("[%s] Dropping rows all-NA in any dataset...", chrom)
  merged <- drop_all_na_in_any_dataset(merged, "dataset", "beta")
  logf("[%s] Post-NA-drop: %d CpGs", chrom, nrow(merged))

  logf("[%s] Imputation (per-CpG, per-condition median):", chrom)
  imp_beta <- impute_median_per_condition(merged, chrom)
  logf("[%s] NAs remaining in imputed beta: %d", chrom, sum(is.na(imp_beta)))

  logf("[%s] Saving to %s...", chrom, out_file)
  saveRDS(
    list(
      chrom    = chrom,
      row_keys = rownames(merged),
      colData  = as.data.frame(colData(merged)),
      beta_imp = imp_beta
    ),
    file = out_file
  )
  rm(merged, imp_beta); invisible(gc(verbose = FALSE))
  logf("[%s] Done.", chrom)
}

# ---- Diagnostics: CpG overlap across conditions (venn, NO MERGE) -----------

diagnose_cpg_overlap_by_condition_nomerge <- function(paths = PATHS$local,
                                                       cfg = FILTER_CFG,
                                                       chrom_list = c(paste0("chr", c(1:22, "X", "Y")))) {

  log_step("CpG overlap diagnostics by condition (no merge)")

  cpg_keys_by_condition <- list(Normal = character(0), Tumour = character(0), metastasis = character(0))

  ds <- load_datasets(paths)

  for (chrom in chrom_list) {
    logf("[%s] Computing per-condition coverage masks (no merge)...", chrom)

    cpgea_chr <- subset_to_chromosome(ds$CPGEA, chrom)
    mcrpc_chr <- subset_to_chromosome(ds$MCRPC, chrom)

    cpgea_chr <- infer_condition(cpgea_chr)
    mcrpc_chr <- infer_condition(mcrpc_chr)

    cov_cpgea <- assay(cpgea_chr, "cov")
    cov_mcrpc <- assay(mcrpc_chr, "cov")

    cond_cpgea <- colData(cpgea_chr)$condition
    cond_mcrpc <- colData(mcrpc_chr)$condition

    keys_cpgea <- coord_key(cpgea_chr)
    keys_mcrpc <- coord_key(mcrpc_chr)

    for (cnd in names(cfg)) {

      if (cnd %in% c("Normal", "Tumour")) {
        cols <- which(cond_cpgea == cnd)
        if (!length(cols)) next
        frac_ok <- DelayedMatrixStats::rowMeans2(
          cov_cpgea[, cols, drop = FALSE] >= cfg[[cnd]]$min_cov, na.rm = TRUE
        )
        this_pass <- !is.nan(frac_ok) & (frac_ok >= cfg[[cnd]]$min_frac)
        cpg_keys_by_condition[[cnd]] <- c(cpg_keys_by_condition[[cnd]], keys_cpgea[this_pass])

      } else if (cnd == "metastasis") {
        cols <- which(cond_mcrpc == cnd)
        if (!length(cols)) next
        frac_ok <- DelayedMatrixStats::rowMeans2(
          cov_mcrpc[, cols, drop = FALSE] >= cfg[[cnd]]$min_cov, na.rm = TRUE
        )
        this_pass <- !is.nan(frac_ok) & (frac_ok >= cfg[[cnd]]$min_frac)
        cpg_keys_by_condition[[cnd]] <- c(cpg_keys_by_condition[[cnd]], keys_mcrpc[this_pass])
      }
    }

    rm(cpgea_chr, mcrpc_chr); invisible(gc(verbose = FALSE))
  }

  rm(ds); invisible(gc(verbose = FALSE))

  counts <- sapply(cpg_keys_by_condition, length)
  logf("CpGs passing coverage per condition - Normal: %d, Tumour: %d, metastasis: %d",
       counts["Normal"], counts["Tumour"], counts["metastasis"])

  common_3 <- Reduce(intersect, cpg_keys_by_condition)
  logf("CpGs common to all 3 conditions: %d", length(common_3))

  all_keys <- unique(unlist(cpg_keys_by_condition))

in_normal     <- all_keys %in% cpg_keys_by_condition$Normal
in_tumour     <- all_keys %in% cpg_keys_by_condition$Tumour
in_metastasis <- all_keys %in% cpg_keys_by_condition$metastasis

membership_count <- in_normal + in_tumour + in_metastasis
at_least_2 <- sum(membership_count >= 2)

logf("CpGs present in at least 2 of 3 conditions: %d", at_least_2)

  if (!requireNamespace("ggVennDiagram", quietly = TRUE)) {
    install.packages("ggVennDiagram", repos = "https://cloud.r-project.org")
  }

  venn_plot <- ggVennDiagram::ggVennDiagram(cpg_keys_by_condition) +
    ggplot2::scale_fill_gradient(low = "white", high = "steelblue") +
    ggplot2::labs(title = "CpGs passing coverage filter by condition (no merge)")

  out_dir <- "results/plots/filtering"
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(file.path(out_dir, "venn_NTM_cpgs_nomerge.pdf"), venn_plot, width = 6, height = 6)

  logf("Venn diagram saved to %s", file.path(out_dir, "venn_NTM_cpgs_nomerge.pdf"))

  invisible(cpg_keys_by_condition)
}

# ---- Reassembly: chromosome-by-chromosome HDF5 write -----------------------

keys_to_granges <- function(keys) {
  m <- do.call(rbind, strsplit(keys, "_"))
  GRanges(seqnames = m[, 1],
          ranges   = IRanges(start = as.integer(m[, 2]),
                             end   = as.integer(m[, 3])))
}

reassemble_to_hdf5 <- function(dir       = CHECKPOINT_DIR,
                               out_h5    = OUT_H5,
                               out_rr    = OUT_ROWRANGES,
                               out_cd    = OUT_COLDATA,
                               chr_order = c(paste0("chr", c(1:22, "X", "Y")))) {

  files <- file.path(dir, sprintf("chr_%s.rds", chr_order))
  missing <- !file.exists(files)
  if (any(missing))
    stop("Missing checkpoints: ", paste(chr_order[missing], collapse = ", "))

  # --- first pass: total row count + colData (cheap) ---
  logf("First pass: counting CpGs and reading colData...")
  obj1     <- readRDS(files[1])
  n_cols   <- nrow(obj1$colData)
  col_data <- obj1$colData
  rm(obj1); invisible(gc(verbose = FALSE))

  n_rows <- sum(sapply(files, function(f) length(readRDS(f)$row_keys)))
  logf("Total: %d CpGs x %d samples (~%.1f GB per assay)",
       n_rows, n_cols, n_rows * n_cols * 8 / 1e9)

  # --- create HDF5 file ---
  if (file.exists(out_h5)) {
    logf("Removing existing %s", out_h5)
    file.remove(out_h5)
  }
  rhdf5::h5createFile(out_h5)
  rhdf5::h5createDataset(out_h5, "beta", dims = c(n_rows, n_cols),
                         chunk = c(1000, n_cols), level = 0,
                         storage.mode = "double")
  rhdf5::h5createDataset(out_h5, "M", dims = c(n_rows, n_cols),
                         chunk = c(1000, n_cols), level = 0,
                         storage.mode = "double")
  logf("HDF5 file created: %s", out_h5)

  # --- second pass: write chromosome by chromosome ---
  row_offset <- 0L
  all_keys   <- character(0)

  for (f in files) {
    chrom <- sub("chr_(chr.+)\\.rds", "\\1", basename(f))
    logf("[%s] Writing to HDF5 (rows %d – %d)...",
         chrom, row_offset + 1L, row_offset + length(readRDS(f)$row_keys))
    obj <- readRDS(f)
    nr  <- nrow(obj$beta_imp)
    idx <- seq(row_offset + 1L, row_offset + nr)

    rhdf5::h5write(obj$beta_imp,            out_h5, "beta", index = list(idx, NULL))
    rhdf5::h5write(beta_to_m(obj$beta_imp), out_h5, "M",    index = list(idx, NULL))

    all_keys   <- c(all_keys, obj$row_keys)
    row_offset <- row_offset + nr
    rm(obj); invisible(gc(verbose = FALSE))
    logf("[%s] Done (offset now %d)", chrom, row_offset)
  }

  rhdf5::H5close()
  logf("HDF5 write complete.")

  # --- save metadata ---
  saveRDS(keys_to_granges(all_keys), out_rr)
  saveRDS(col_data, out_cd)
  logf("Metadata saved: %s, %s", out_rr, out_cd)
}

#' Load the reassembled SE from HDF5 + metadata RDS files.
load_integrated_se <- function(out_h5   = OUT_H5,
                               out_rr   = OUT_ROWRANGES,
                               out_cd   = OUT_COLDATA) {
  beta     <- HDF5Array::HDF5Array(out_h5, "beta")
  M        <- HDF5Array::HDF5Array(out_h5, "M")
  col_data <- readRDS(out_cd)
  row_gr   <- readRDS(out_rr)
  SummarizedExperiment(
    assays    = list(beta = beta, M = M),
    rowRanges = row_gr,
    colData   = S4Vectors::DataFrame(col_data)
  )
}

# ---- Driver -----------------------------------------------------------------

run_pipeline <- function(paths = PATHS$hpc, only_chrom = NULL) {
  log_start()
  ensure_dirs()

  if (identical(only_chrom, "REASSEMBLE")) {
    log_step("Reassemble to HDF5")
    reassemble_to_hdf5()
    logf("Done. Total wall time: %.1f min",
         as.numeric(difftime(Sys.time(), .PIPE_T0, units = "mins")))
    return(invisible(NULL))
  }

  if (identical(only_chrom, "DIAGNOSE_VENN_NOMERGE")) {
    log_step("CpG overlap diagnostics (venn, no merge)")
    diagnose_cpg_overlap_by_condition_nomerge(paths = paths)
    logf("Done. Total wall time: %.1f min",
         as.numeric(difftime(Sys.time(), .PIPE_T0, units = "mins")))
    return(invisible(NULL))
  }

  chrom_list <- if (!is.null(only_chrom)) only_chrom
                else c(paste0("chr", c(2, 1, 3:22, "X", "Y")))

  for (chrom in chrom_list) {
    process_chromosome(chrom, paths, CHECKPOINT_DIR, FILTER_CFG)
  }

  if (is.null(only_chrom)) {
    log_step("Reassemble to HDF5")
    reassemble_to_hdf5()
    logf("Done. Total wall time: %.1f min",
         as.numeric(difftime(Sys.time(), .PIPE_T0, units = "mins")))
    return(invisible(NULL))
  }

  logf("Finished chromosome '%s'. Total wall time: %.1f min",
       only_chrom, as.numeric(difftime(Sys.time(), .PIPE_T0, units = "mins")))
  invisible(NULL)
}

# run_pipeline(PATHS$local)
run_pipeline(PATHS$local, only_chrom = "DIAGNOSE_VENN_NOMERGE")


# ---- Entry point ------------------------------------------------------------
# All chromosomes, one session:
# run_pipeline(PATHS$local)
#   se <- load_integrated_se()
#
# SLURM array (one chrom per task):
#   args  <- commandArgs(trailingOnly = TRUE)
#   chrom <- if (length(args)) args[[1]] else NULL
#   run_pipeline(PATHS$hpc, only_chrom = chrom)
#
# After all array tasks finish:
#   run_pipeline(PATHS$hpc, only_chrom = "REASSEMBLE")
#   se <- load_integrated_se()
#
# Diagnostics venn (CpG overlap by condition, no merge, all chromosomes):
#   run_pipeline(PATHS$local, only_chrom = "DIAGNOSE_VENN_NOMERGE")