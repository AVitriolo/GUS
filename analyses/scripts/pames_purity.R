suppressPackageStartupMessages({
  library(HDF5Array)
  library(SummarizedExperiment)
  library(GenomicRanges)
  library(DelayedMatrixStats)
  library(PAMES)
})

options(scipen = 999)

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

cpgea_path      <- args$cpgea
min_cov         <- as.numeric(args$min_cov)
min_sample_frac <- as.numeric(args$min_sample_frac)
min_mean_diff   <- as.numeric(args$min_mean_diff)
max_na_data     <- as.numeric(args$max_na_data)
percentile_low  <- as.numeric(args$percentile_low)
percentile_high <- as.numeric(args$percentile_high)
max_sites       <- as.integer(args$max_sites)
min_distance    <- as.numeric(args$min_distance)
cores           <- as.integer(args$cores)
out_purity      <- args$output_purity
out_sites       <- args$output_info_sites

se <- HDF5Array::loadHDF5SummarizedExperiment(cpgea_path)

condition <- SummarizedExperiment::colData(se)$condition
assertthat::assert_that(all(condition %in% c("Normal", "Tumour")),
                        msg = "Expected only Normal/Tumour samples; metastatic samples must not be passed in.")

is_tumor   <- condition == "Tumour"
is_control <- condition == "Normal"

message(sprintf("[%s] Samples: %d Tumour, %d Normal", Sys.time(), sum(is_tumor), sum(is_control)))

## Two-stage candidate CpG filter, computed blockwise (HDF5-backed) so the
## full 29M-CpG x 374-sample matrix is never materialised in memory. NAs are
## handled via na.rm in the blockwise reductions rather than by writing zeros
## into the DelayedMatrix, which would force an expensive full materialisation.

## Stage 1: coverage. Require >= min_cov reads in >= min_sample_frac of samples,
## in both groups independently.
cov <- SummarizedExperiment::assay(se, "cov")

message(sprintf("[%s] Filtering CpGs by coverage (>= %d reads in >= %.0f%% of samples, in both groups)...",
                Sys.time(), min_cov, 100 * min_sample_frac))

tumor_covered   <- DelayedMatrixStats::rowSums2(cov[, is_tumor,   drop = FALSE] >= min_cov, na.rm = TRUE)
control_covered <- DelayedMatrixStats::rowSums2(cov[, is_control, drop = FALSE] >= min_cov, na.rm = TRUE)

keep_cov <- which(tumor_covered   >= min_sample_frac * sum(is_tumor) &
                  control_covered >= min_sample_frac * sum(is_control))

message(sprintf("[%s] Candidate CpGs after coverage filter: %d / %d",
                Sys.time(), length(keep_cov), nrow(se)))
assertthat::assert_that(length(keep_cov) > 0, msg = "No CpG passed the coverage filter.")

rm(cov, tumor_covered, control_covered)
gc()

## Stage 2: mean methylation difference. Still blockwise/HDF5-backed, applied
## only to the coverage survivors, to discard CpGs with no tumor/normal signal
## before the informative-site search (which would reject them anyway).
se_cov <- se[keep_cov, ]
beta_cov <- SummarizedExperiment::assay(se_cov, "beta")

message(sprintf("[%s] Filtering CpGs by mean methylation difference (>= %.2f)...",
                Sys.time(), min_mean_diff))

tumor_mean   <- DelayedMatrixStats::rowMeans2(beta_cov[, is_tumor,   drop = FALSE], na.rm = TRUE)
control_mean <- DelayedMatrixStats::rowMeans2(beta_cov[, is_control, drop = FALSE], na.rm = TRUE)

keep_diff <- which(abs(tumor_mean - control_mean) >= min_mean_diff)

message(sprintf("[%s] Candidate CpGs after mean-difference filter: %d / %d",
                Sys.time(), length(keep_diff), length(keep_cov)))
assertthat::assert_that(length(keep_diff) > 0, msg = "No CpG passed the mean-difference filter.")

rm(beta_cov, tumor_mean, control_mean)
gc()

## Realise only the final filtered subset into memory
sub  <- se_cov[keep_diff, ]
rm(se, se_cov)
gc()
beta <- as.matrix(SummarizedExperiment::assay(sub, "beta"))
rownames(beta) <- paste0(GenomicRanges::seqnames(SummarizedExperiment::rowRanges(sub)), "_",
                         GenomicRanges::start(SummarizedExperiment::rowRanges(sub)))

tumor_table   <- beta[, is_tumor,   drop = FALSE]
control_table <- beta[, is_control, drop = FALSE]

ref_table <- data.frame(
  Chromosome         = as.character(GenomicRanges::seqnames(SummarizedExperiment::rowRanges(sub))),
  Genomic_Coordinate = GenomicRanges::start(SummarizedExperiment::rowRanges(sub))
)

rm(sub)
gc()

## PAMES site-based workflow: AUC -> informative sites -> purity
## NOTE: PAMES::get_AUC has an inverted comparison (keeps sites with NA-fraction
## >= max_NA_data instead of <=), so any non-zero max_na_data silently drops every
## site when data is fully complete (as ours is, by construction of the coverage
## prefilter above). max_na_data must stay 0 for that filter to be a no-op instead
## of discarding everything.
auc <- get_AUC(tumor_table, control_table, max_NA_data = max_na_data, cores = cores)

## PAMES::find_informative_sites also has a bug: the line meant to put hypo-site
## AUC on a comparable "higher = stronger separation" scale as hyper sites
## (`diff_meth_sites_df$AUC[diff_meth_sites_df$Site_type == "Hypo"] <- 1-AUC`)
## checks for "Hypo" (capital H) but Site_type values are lowercase "hypo", so it
## never fires. Sites are then ranked by that un-flipped, raw AUC, which means
## selection (and the min_distance dedup, which consumes rows in ranked order)
## picks the *weakest* hypo sites (AUC near 0.2) before the strongest (AUC near
## 0), the opposite of the intended "most informative first". Hyper-site ranking
## is unaffected (their AUC is already "higher = stronger", untouched by the bug).
##
## Fix: get every candidate via full_info=TRUE (skips PAMES's own buggy
## selection), re-rank hypo sites by (1 - AUC) ourselves, then reuse PAMES's own
## (correct) distance-dedup helper on our correctly-sorted candidates.
all_sites <- find_informative_sites(tumor_table, control_table, auc, ref_table,
                                    cores = cores, max_sites = max_sites,
                                    min_distance = min_distance,
                                    percentiles = c(percentile_low, percentile_high),
                                    full_info = TRUE)

select_ranked <- function(df, score, n) {
  if (nrow(df) == 0) return(df[0, ])
  df <- df[order(score, decreasing = TRUE), ]
  df$Keep_site <- PAMES:::remove_close_sites(df, n, min_distance)
  head(df[df$Keep_site, ], n)
}

hyper_df <- all_sites[all_sites$Site_type == "hyper", ]
hypo_df  <- all_sites[all_sites$Site_type == "hypo",  ]

hyper_kept <- select_ranked(hyper_df, hyper_df$AUC,         max_sites / 2)
hypo_kept  <- select_ranked(hypo_df,  1 - hypo_df$AUC,      max_sites / 2)

info_sites <- rbind(hyper_kept, hypo_kept)[, c("Probe", "Site_type")]

message(sprintf("[%s] Selected %d hyper- and %d hypo-methylated sites (bug-corrected ranking)",
                Sys.time(), nrow(hyper_kept), nrow(hypo_kept)))

assertthat::assert_that(nrow(info_sites) > 0,
                        msg = "PAMES did not find any informative site; consider relaxing filters.")

full_table <- cbind(tumor_table, control_table)
purity_df <- get_purity(full_table, info_sites)
purity_df$Condition <- ifelse(purity_df$Sample %in% colnames(tumor_table), "Tumour", "Normal")

dir.create(dirname(out_purity), recursive = TRUE, showWarnings = FALSE)
write.table(purity_df, file = out_purity, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(info_sites, file = out_sites, sep = "\t", quote = FALSE, row.names = FALSE)

message(sprintf("[%s] Done. Purity written to %s (%d sites written to %s)",
                Sys.time(), out_purity, nrow(info_sites), out_sites))
