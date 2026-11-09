# =============================================================================
# inhouse_delta_methylation.R
# Metastasis vs tumour deltas of the in-house WGBS cohort (methrix HDF5 SE,
# beta + cov, no imputation, no ComBat) on the methyldriver regions
# (md_percpg BED: chr start end CpGID clusterID; CpGID = chr_pos).
#   values : beta where cov >= --min_cov, else NA
#   CpG kept if >= --min_frac of the (non-excluded) pairs have BOTH samples
#            covered >= --min_cov; region kept if >= 1 of its CpGs is kept
#   per kept CpG, over the kept CpGs of each region (mean of the CpG values):
#     paired   : delta = mean over complete pairs of (M_i - T_i)
#                t     = mean(d) / (sd(d) / sqrt(n_pairs)), CpG-level, then
#                        averaged per region
#     unpaired : delta = mean(M) - mean(T), all samples with cov >= min_cov
# Rows are matched to the SE by chr + pos (SE row order is not the BED order).
# Output (prefix): _paired.tsv, _paired_t.tsv, _unpaired.tsv (clusterID chr
# start end n_CpGs n_CpGs_kept T M delta [n_pairs]) and _regions.tsv (the
# md_regions rows of the kept regions, for methyldriver_w20.R).
# =============================================================================
suppressPackageStartupMessages({
  library(SummarizedExperiment); library(HDF5Array); library(data.table); library(matrixStats)
})
args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
MIN_COV  <- as.integer(if (is.null(args$min_cov)) 5 else args$min_cov)
MIN_FRAC <- as.numeric(if (is.null(args$min_frac)) 0.8 else args$min_frac)
exclude  <- if (is.null(args$exclude_pairs)) character(0) else strsplit(args$exclude_pairs, ",")[[1]]
prefix   <- args$out_prefix
dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)

pairs <- fread(args$pairs)                       # pair, tumour, metastasis
pairs <- pairs[!pair %in% exclude]
cat(sprintf("pairs used: %d (%s)\n", nrow(pairs), paste(pairs$pair, collapse = ",")))

se <- loadHDF5SummarizedExperiment(args$se_dir)
stopifnot(all(c(pairs$tumour, pairs$metastasis) %in% colnames(se)))
bed <- fread(args$bed, header = FALSE, col.names = c("chr", "start", "end", "CpGID", "clusterID"))
bed[, pos := as.integer(sub(".*_", "", CpGID))]
cpgs <- unique(bed[, .(CpGID, chr, pos)])
cat(sprintf("BED: %d CpG rows, %d unique CpGs, %d regions\n", nrow(bed), nrow(cpgs), uniqueN(bed$clusterID)))

chr_all <- as.character(seqnames(se)); pos_all <- start(se)
res <- list()
for (ch in unique(cpgs$chr)) {
  k <- cpgs[chr == ch]
  idx <- which(chr_all == ch)
  r <- idx[match(k$pos, pos_all[idx])]
  if (anyNA(r)) stop(sprintf("%s: %d BED CpGs not in the SE", ch, sum(is.na(r))))
  o  <- order(r)
  cv <- as.matrix(assay(se, "cov")[r[o], c(pairs$tumour, pairs$metastasis), drop = FALSE])
  b  <- as.matrix(assay(se, "beta")[r[o], c(pairs$tumour, pairs$metastasis), drop = FALSE])
  b[is.na(cv) | cv < MIN_COV] <- NA
  np <- nrow(pairs)
  bt <- b[, seq_len(np), drop = FALSE]; bm <- b[, np + seq_len(np), drop = FALSE]
  d  <- bm - bt                                                        # NA unless both covered
  n_pairs <- rowSums(!is.na(d))
  mean_d  <- rowMeans(d, na.rm = TRUE)
  sd_d    <- rowSds(d, na.rm = TRUE)
  res[[ch]] <- data.table(CpGID = k$CpGID[o], keep = n_pairs / np >= MIN_FRAC, n_pairs = n_pairs,
                          paired = mean_d, paired_t = mean_d / (sd_d / sqrt(n_pairs)),
                          T = rowMeans(bt, na.rm = TRUE), M = rowMeans(bm, na.rm = TRUE))
}
cpg <- rbindlist(res)
cpg[, unpaired := M - T]
# a CpG with identical differences in all pairs has sd 0: t is undefined, keep it out of the t mean
cpg[!is.finite(paired_t), paired_t := NA]
cat(sprintf("CpGs kept (>= %.0f%% of %d pairs with both samples >= %dx): %d of %d\n",
            100 * MIN_FRAC, nrow(pairs), MIN_COV, sum(cpg$keep), nrow(cpg)))

x <- merge(bed, cpg, by = "CpGID")[keep == TRUE]
reg <- x[, .(chr = chr[1], start = min(start), end = max(end), n_CpGs_kept = .N,
             T = mean(T), M = mean(M), n_pairs = round(mean(n_pairs), 1),
             paired = mean(paired), paired_t = mean(paired_t, na.rm = TRUE), unpaired = mean(unpaired)),
         by = clusterID]
reg <- merge(reg, bed[, .(n_CpGs = .N), by = clusterID], by = "clusterID")
reg <- reg[is.finite(paired_t)]
cat(sprintf("regions kept: %d of %d\n", nrow(reg), uniqueN(bed$clusterID)))

for (v in c("paired", "paired_t", "unpaired")) {
  out <- reg[, .(clusterID, chr, start, end, n_CpGs, n_CpGs_kept, T, M, delta = get(v), n_pairs)][order(-abs(delta))]
  fwrite(out, sprintf("%s_%s.tsv", prefix, v), sep = "\t")
}
md <- fread(args$md_regions)
fwrite(md[mdID %in% reg$clusterID], sprintf("%s_regions.tsv", prefix), sep = "\t")
cat(sprintf("cor(paired, unpaired) over regions: %.3f; median |paired - unpaired|: %.4f\n",
            cor(reg$paired, reg$unpaired), median(abs(reg$paired - reg$unpaired))))
cat("written", prefix, "_{paired,paired_t,unpaired,regions}.tsv\n")
