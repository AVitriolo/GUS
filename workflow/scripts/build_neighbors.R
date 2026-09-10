#!/usr/bin/env Rscript
set.seed(123); options(scipen = 999)
source("workflow/scripts/helpers/get_encode_info.R")
source("workflow/scripts/helpers/make_feature_matrix_from_bigwigs_and_bed.R")
source("workflow/scripts/helpers/most_similar_n_regions.R")

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
input_bed        <- args$input_bed
bigwig_dir       <- args$bigwig_dir
clusters_to_keep <- args$clusters_to_keep
output_neighbors <- args$output_neighbors     # .rds
output_matrix    <- args$output_matrix        # .rds
num_cores        <- if (!is.null(args$num_cores)) as.integer(args$num_cores) else 8
n_neighbors      <- if (!is.null(args$n)) as.integer(args$n) else 100
n_pcs            <- if (!is.null(args$n_pcs)) as.integer(args$n_pcs) else 20
exclude_files    <- if (!is.null(args$exclude_files)) strsplit(args$exclude_files, ",")[[1]] else character(0)

for (f in c(output_neighbors, output_matrix))
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)

message("== building feature matrix ==")
fm <- make_feature_matrix_from_bigwigs_and_bed(
  bigwig_dir = bigwig_dir, genomic_regions_bed = input_bed,
  exclude_files = exclude_files, num_cores = num_cores)

if (!is.null(clusters_to_keep)) {
  keep <- read.table(clusters_to_keep, header = FALSE, stringsAsFactors = FALSE)$V1
  fm <- fm[rownames(fm) %in% keep, , drop = FALSE]
}
message(sprintf("feature matrix: %d regions x %d features", nrow(fm), ncol(fm)))

message("== finding most-similar regions ==")
nb <- most_similar_n_regions(fm, n = n_neighbors, n_pcs = n_pcs)

saveRDS(fm, output_matrix)
saveRDS(nb, output_neighbors)
message("done")