options(scipen = 999)

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

input_path_corr     <- args$input_path_corr
probe_bed_path       <- args$probe_bed
TxID                 <- args$tx_id
output_path_cpg_list <- args$output_path_cpg_list
output_path_corr     <- args$output_path_corr

# ---- Illumina HM450 probe CpGIDs (chr_start, matching the pipeline's CpGID convention) ----
# The probe BED is 0-based half-open (standard BED convention), but the pipeline's
# CpGID / rowRanges are 1-based (confirmed against integrated_wgbs_combat_rowRanges.rds:
# BED row "chr16 53434199 53434201" for probe cg00000029 corresponds to the 1-based
# CpG position 53434200, i.e. BED start + 1) - so the +1 below is required, not cosmetic.
probes <- read.table(probe_bed_path, header = FALSE, sep = "\t",
                      col.names = c("chr", "start", "end", "name"))
probe_ids <- paste0(probes$chr, "_", probes$start + 1)

# ---- CpGs already correlated for this TxID (data/corr_NTM/..._{TxID}) ----
corr <- read.delim(input_path_corr, check.names = FALSE, header = TRUE, sep = "\t")
all_ids  <- colnames(corr)
keep_ids <- intersect(all_ids, probe_ids)

cat(sprintf("[%s] %d CpGs -> %d Illumina HM450 probe CpGs\n", TxID, length(all_ids), length(keep_ids)))

if (length(keep_ids) < 2) {
  cat(sprintf("[%s] only %d probe CpG(s) survive - skipping (need >= 2 for clustering).\n",
              TxID, length(keep_ids)))
  quit(save = "no", status = 1)
}

# ---- CpG list for unpack_RSE_NTM.R (same 2-col TxID<TAB>CpGID format as unpack_topK_NTM) ----
dir.create(dirname(output_path_cpg_list), recursive = TRUE, showWarnings = FALSE)
writeLines(paste(TxID, keep_ids, sep = "\t"), output_path_cpg_list)

# ---- correlation matrix subset to the surviving CpGs (pairwise correlations are unaffected
#      by dropping other CpGs, so this is equivalent to recomputing from scratch) ----
corr_filtered <- corr[keep_ids, keep_ids]

dir.create(dirname(output_path_corr), recursive = TRUE, showWarnings = FALSE)
write.table(corr_filtered, file = output_path_corr, sep = "\t",
            quote = FALSE, row.names = TRUE, col.names = TRUE)

cat(sprintf("[%s] written %d x %d probe-restricted correlation matrix\n",
            TxID, nrow(corr_filtered), ncol(corr_filtered)))
