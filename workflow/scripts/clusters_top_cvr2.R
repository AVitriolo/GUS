#!/usr/bin/env Rscript
# =============================================================================
# clusters_from_top_cvr2.R
# Build the cluster-ID keep-list for downstream (MethylDriver / gimme):
# clusters whose TRANSCRIPT has cv_r2 >= the 3rd quartile of the cv_r2
# distribution across all modeled transcripts.
#
# Usage:
#   Rscript clusters_from_top_cvr2.R \
#     --perf_dir=results/performance \
#     --bed=data/to_methyldriver/selected_features_hg38_v38_pc_50000_25_10_3.bed \
#     --out=resources/clusters_top_cvr2.txt
# =============================================================================
args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
perf_dir <- args$perf_dir
bed_path <- args$bed
out_path <- args$out
quant    <- if (!is.null(args$quantile)) as.numeric(args$quantile) else 0.75

# ---- read cv_r2 per transcript from the performance files ----
files <- list.files(perf_dir, full.names = TRUE)
files <- files[!file.info(files)$isdir]
cat(sprintf("performance files: %d\n", length(files)))

read_cvr2 <- function(f) {
    df <- tryCatch(read.table(f, sep = "\t", header = FALSE,
                              stringsAsFactors = FALSE, row.names = 1),
                   error = function(e) NULL)
    if (is.null(df) || !("cv_r2" %in% rownames(df))) return(NULL)
    txid <- sub(".*_(ENST[0-9]+)$", "\\1", basename(f))
    data.frame(TxID = txid, cv_r2 = as.numeric(df["cv_r2", 1]),
               stringsAsFactors = FALSE)
}
perf <- do.call(rbind, lapply(files, read_cvr2))
perf <- perf[!is.na(perf$cv_r2), ]
cat(sprintf("transcripts with cv_r2: %d\n", nrow(perf)))

# ---- threshold = 3rd quartile of the cv_r2 distribution ----
thr <- quantile(perf$cv_r2, probs = quant, na.rm = TRUE)
top_tx <- perf$TxID[perf$cv_r2 >= thr]
cat(sprintf("cv_r2 %.0f%% quantile = %.4f ; transcripts kept: %d of %d\n",
            100 * quant, thr, length(top_tx), nrow(perf)))

# ---- map to cluster IDs: keep BED clusters whose TxID is in top_tx ----
bed <- read.table(bed_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
cluster_ids <- bed$V4                              # <TxID>_<cluster_id>
cluster_tx  <- sub("_[^_]+$", "", cluster_ids)     # strip trailing _<cluster_id>
keep <- cluster_ids[cluster_tx %in% top_tx]
keep <- unique(keep)

cat(sprintf("clusters kept: %d of %d\n", length(keep), length(unique(cluster_ids))))

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
writeLines(keep, out_path)
cat(sprintf("written -> %s\n", out_path))