options(scipen = 999)

# For every existing cluster (CpGID <TAB> cluster_id file produced by merge_CpGs),
# measure how much of it is made up of Illumina HM450 probe CpGs: per-cluster and
# per-TxID probe coverage. Read-only diagnostic - does not modify or regenerate
# any clustering output, and works directly against whichever clusters_maps
# directory you point it at (the existing full-CpG clustering by default, or the
# HM450-restricted one from run_hm450_validation.sh, where coverage is 100% by
# construction and mainly useful as a sanity check).
#
# Reads all cluster map files first, then does a single vectorized probe lookup
# (rather than one match() per file) so the ~486k-probe hash table is built once,
# not once per TxID - matters at the ~27k-TxID scale this pipeline runs at.
#
# Usage:
#   Rscript analyses/scripts/cluster_probe_coherence.R \
#     --clusters_dir=data/clustering/clusters_maps \
#     --probe_bed=resources/illumina_probes/HM450_hg38_probes.bed \
#     --output_path_clusters=results/cluster_probe_coherence_clusters.tsv \
#     --output_path_summary=results/cluster_probe_coherence_summary.tsv

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

clusters_dir          <- args$clusters_dir
probe_bed_path        <- args$probe_bed
output_path_clusters  <- args$output_path_clusters
output_path_summary   <- args$output_path_summary

# ---- Illumina HM450 probe CpGIDs ----
# The probe BED is 0-based half-open (standard BED convention); the pipeline's
# CpGID / rowRanges are 1-based, hence the +1 (verified against
# integrated_wgbs_combat_rowRanges.rds: BED row "chr16 53434199 53434201" for
# probe cg00000029 is the 1-based CpG position 53434200 = BED start + 1).
probes <- data.table::fread(probe_bed_path, header = FALSE,
                             col.names = c("chr", "start", "end", "name"))
probe_ids <- probes[, paste0(chr, "_", start + 1)]
cat(sprintf("Loaded %d Illumina HM450 probe positions\n", length(probe_ids)))

files <- list.files(clusters_dir, pattern = "_clusters_map\\.txt$", full.names = TRUE)
if (length(files) == 0)
  stop(sprintf("No *_clusters_map.txt files found in %s", clusters_dir))
cat(sprintf("Found %d cluster map files in %s - reading...\n", length(files), clusters_dir))

read_one <- function(f) {
  TxID <- sub("_clusters_map\\.txt$", "", basename(f))
  TxID <- sub(".*_", "", TxID)
  dt <- tryCatch(
    data.table::fread(f, header = FALSE, col.names = c("CpGID", "cluster_id")),
    error = function(e) NULL
  )
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  dt[, TxID := TxID]
  dt
}

chunks <- vector("list", length(files))
for (i in seq_along(files)) {
  chunks[[i]] <- read_one(files[i])
  if (i %% 2000 == 0) cat(sprintf("  read %d / %d\n", i, length(files)))
}
all_cpgs <- data.table::rbindlist(chunks, fill = TRUE)
rm(chunks)
cat(sprintf("Combined %d CpG-cluster rows across %d TxIDs\n",
            nrow(all_cpgs), length(unique(all_cpgs$TxID))))

# single vectorized lookup - builds the probe hash table once
all_cpgs[, is_probe := CpGID %in% probe_ids]

cluster_out <- all_cpgs[, .(n_CpGs = .N, n_probe_CpGs = sum(is_probe)), by = .(TxID, cluster_id)]
cluster_out[, frac_probe_CpGs := round(n_probe_CpGs / n_CpGs, 4)]

summary_out <- cluster_out[, .(
  n_CpGs                = sum(n_CpGs),
  n_probe_CpGs           = sum(n_probe_CpGs),
  n_clusters             = .N,
  n_clusters_with_probe  = sum(n_probe_CpGs > 0),
  n_clusters_all_probe   = sum(n_probe_CpGs == n_CpGs)
), by = TxID]
summary_out[, frac_probe_CpGs := round(n_probe_CpGs / n_CpGs, 4)]

dir.create(dirname(output_path_clusters), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_path_summary),  recursive = TRUE, showWarnings = FALSE)
data.table::fwrite(cluster_out, output_path_clusters, sep = "\t")
data.table::fwrite(summary_out, output_path_summary,  sep = "\t")

cat(sprintf(
  "\n%d TxIDs, %d clusters total\n", nrow(summary_out), sum(summary_out$n_clusters)))
cat(sprintf(
  "Overall: %.2f%% of all clustered CpGs are Illumina HM450 probes (%d / %d)\n",
  100 * sum(summary_out$n_probe_CpGs) / sum(summary_out$n_CpGs),
  sum(summary_out$n_probe_CpGs), sum(summary_out$n_CpGs)))
cat(sprintf(
  "Median per-TxID probe coverage: %.2f%%\n", 100 * stats::median(summary_out$frac_probe_CpGs)))
cat(sprintf(
  "Clusters containing >=1 probe CpG: %.2f%% (%d / %d)\n",
  100 * sum(summary_out$n_clusters_with_probe) / sum(summary_out$n_clusters),
  sum(summary_out$n_clusters_with_probe), sum(summary_out$n_clusters)))
cat(sprintf(
  "Clusters made ENTIRELY of probe CpGs: %.2f%% (%d / %d)\n",
  100 * sum(summary_out$n_clusters_all_probe) / sum(summary_out$n_clusters),
  sum(summary_out$n_clusters_all_probe), sum(summary_out$n_clusters)))

cat(sprintf("\nWritten: %s (per-cluster)\n", output_path_clusters))
cat(sprintf("Written: %s (per-TxID summary)\n", output_path_summary))
