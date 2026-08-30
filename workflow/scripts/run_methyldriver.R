#!/usr/bin/env Rscript
# =============================================================================
# run_methyldriver.R
# End-to-end MethylDriver stage for the selected clusters:
#   1. build the region x chromatin-feature matrix from ENCODE bigWigs over the
#      selected-features BED  (make_feature_matrix_from_bigwigs_and_bed)
#   2. find each cluster's most-similar neighbours in feature space
#      (most_similar_n_regions)
#   3. run MethylDriver on the per-cluster methylation change and draw a volcano
#
# --clusters_to_keep is OPTIONAL: if given, restricts to those cluster IDs
# (e.g. clusters with >= 3 CpGs); if omitted, all clusters are used.
# =============================================================================
set.seed(123)
options(scipen = 999)

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(ggrepel)
})

source("workflow/scripts/helpers/get_encode_info.R")
source("workflow/scripts/helpers/make_feature_matrix_from_bigwigs_and_bed.R")
source("workflow/scripts/helpers/most_similar_n_regions.R")
source("workflow/scripts/helpers/MethylDriver.R")

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

input_bed         <- args$input_bed
input_meth_diff   <- args$input_meth_diff
bigwig_dir        <- args$bigwig_dir
conversion_table  <- args$conversion_table
clusters_to_keep  <- args$clusters_to_keep
output_results    <- args$output_results
output_volcano    <- args$output_volcano
output_full_table <- args$output_full_table
num_cores         <- if (!is.null(args$num_cores)) as.integer(args$num_cores) else 8
n_neighbors       <- if (!is.null(args$n)) as.integer(args$n) else 100
n_pcs             <- if (!is.null(args$n_pcs)) as.integer(args$n_pcs) else 20
q_cutoff          <- if (!is.null(args$q)) as.numeric(args$q) else 0.05
bigwig_pattern    <- args$bigwig_pattern
exclude_files     <- if (!is.null(args$exclude_files)) strsplit(args$exclude_files, ",")[[1]] else character(0)

for (f in c(output_results, output_volcano, output_full_table))
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)

# ---- 1. feature matrix -------------------------------------------------------
message("== building feature matrix ==")
feature_matrix_ml_clusters <- make_feature_matrix_from_bigwigs_and_bed(
  bigwig_dir          = bigwig_dir,
  genomic_regions_bed = input_bed,
  bigwig_pattern      = bigwig_pattern,
  exclude_files       = exclude_files,
  num_cores           = num_cores)

# ---- per-cluster methylation change -----------------------------------------
# cluster_delta_methylation.R writes columns: clusterID, ..., delta
mean_meth_diff <- read.table(input_meth_diff, header = TRUE, sep = "\t",
                             stringsAsFactors = FALSE)
home_region_methyl_change <- setNames(mean_meth_diff$delta,
                                      mean_meth_diff$clusterID)

feature_matrix_ml_clusters.f <-
  feature_matrix_ml_clusters[rownames(feature_matrix_ml_clusters) %in%
                               mean_meth_diff$clusterID, , drop = FALSE]

# ---- OPTIONAL cluster filter -------------------------------------------------
if (!is.null(clusters_to_keep)) {
  keep <- read.table(clusters_to_keep, header = FALSE, stringsAsFactors = FALSE)$V1
  message(sprintf("clusters_to_keep: restricting to %d clusters", length(keep)))
  feature_matrix_ml_clusters.ff <-
    feature_matrix_ml_clusters.f[rownames(feature_matrix_ml_clusters.f) %in% keep, , drop = FALSE]
} else {
  message("clusters_to_keep: not given - using all clusters")
  feature_matrix_ml_clusters.ff <- feature_matrix_ml_clusters.f
}
message(sprintf("feature matrix: %d clusters x %d features",
                nrow(feature_matrix_ml_clusters.ff), ncol(feature_matrix_ml_clusters.ff)))

# ---- 2. nearest neighbours in feature space ---------------------------------
message("== finding most-similar regions ==")
closest_regions_list_ML.ff <- most_similar_n_regions(
  feature_matrix_ml_clusters.ff, n = n_neighbors, n_pcs = n_pcs)

# ---- align the three objects to the neighbour-list order --------------------
target_order <- names(closest_regions_list_ML.ff)

home_region_methyl_change.ff <- home_region_methyl_change[target_order]

transcript_ids <- sub("_.*$", "", target_order)
conversion_tbl <- read.table(conversion_table, header = TRUE, sep = "\t",
                             stringsAsFactors = FALSE)
gene_names   <- conversion_tbl$gene_name[match(transcript_ids, conversion_tbl$transcript_id)]
regions_genes <- setNames(gene_names, target_order)

stopifnot(all(names(home_region_methyl_change.ff) == target_order))
stopifnot(all(names(regions_genes) == target_order))

# ---- 3. MethylDriver --------------------------------------------------------
message("== running MethylDriver ==")
results_methyldriver <- MethylDriver(
  home_region_methyl_change = home_region_methyl_change.ff,
  regions_genes             = regions_genes,
  closest_regions_list      = closest_regions_list_ML.ff,
  q_value_cutoff            = q_cutoff)

saveRDS(results_methyldriver, file = output_results)
write.table(results_methyldriver$full_results, file = output_full_table,
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- volcano ----------------------------------------------------------------
df <- results_methyldriver$full_results
df$significance <- ifelse(df$q_value < q_cutoff & df$home_region_methyl_change > 0, "Hyper",
                   ifelse(df$q_value < q_cutoff & df$home_region_methyl_change < 0, "Hypo", "NS"))

top_genes <- df[order(df$p_value), ][seq_len(min(55, nrow(df))), ]

pdf(output_volcano, width = 12, height = 12)
print(
  ggplot(df, aes(x = home_region_methyl_change, y = -log10(p_value), color = significance)) +
    geom_point(alpha = 0.7) +
    ggrepel::geom_text_repel(data = top_genes, aes(label = gene),
                             color = "black", size = 6,
                             box.padding = 0.3, point.padding = 0.2, max.overlaps = 100) +
    scale_color_manual(values = c("Hypo" = "steelblue3", "Hyper" = "#F8766B", "NS" = "grey")) +
    labs(x = "Methylation Change", y = "-log10(p-value)",
         title = "Methylation Changes - MethylDriver", color = "Significance") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 22),
          axis.title = element_text(size = 20), axis.text = element_text(size = 18),
          legend.title = element_blank(), legend.text = element_text(size = 18),
          legend.position = "top")
)
dev.off()

message(sprintf("done: %d hyper, %d hypo genes",
                length(results_methyldriver$hypermethylated_genes),
                length(results_methyldriver$hypomethylated_genes)))