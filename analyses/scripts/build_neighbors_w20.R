#!/usr/bin/env Rscript
# =============================================================================
# build_neighbors_w20.R  (run on the HPC, where the bigWigs are)
# w20 counterpart of build_neighbors.R:
#   1. feature matrix over ALL distinct windows (--input_bed), with
#      make_feature_matrix_from_bigwigs_and_bed() as build_neighbors.R, so a
#      different cv_r2 cutoff later needs no new bigWig pass
#   2. neighbours of the methyldriver regions (--md_regions: models above the
#      mean cv_r2, identical transcript regions collapsed), skipping neighbours
#      closer than --min_dist bp (helpers/neighbors_min_dist.R)
# Outputs are the files run_methyldriver_gimme_w20.sh / methyldriver_w20.R
# reuse: featmatrix_<tag>.rds and neighbors_<tag>.rds.
# Column targets come from --columns (names of the previous matrix,
# <accession>_<target>) instead of the ENCODE portal, so the run does not
# depend on internet access from the compute node.
# If --output_matrix exists it is reused (only the neighbours are redone).
# =============================================================================
set.seed(123); options(scipen = 999)
suppressPackageStartupMessages(library(data.table))
source("workflow/scripts/helpers/make_feature_matrix_from_bigwigs_and_bed.R")
source("analyses/scripts/helpers/neighbors_min_dist.R")

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
input_bed        <- args$input_bed
md_path          <- args$md_regions
bigwig_dir       <- args$bigwig_dir
output_matrix    <- args$output_matrix
output_neighbors <- args$output_neighbors
cols             <- readLines(args$columns)
exclude_files    <- if (!is.null(args$exclude_files)) strsplit(args$exclude_files, ",")[[1]] else character(0)
num_cores        <- as.integer(if (is.null(args$num_cores)) 8 else args$num_cores)
n_neighbors      <- as.integer(if (is.null(args$n)) 100 else args$n)
n_pcs            <- as.integer(if (is.null(args$n_pcs)) 20 else args$n_pcs)
min_dist         <- as.integer(if (is.null(args$min_dist)) 1000 else args$min_dist)
for (f in c(output_matrix, output_neighbors)) dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)

target <- setNames(sub("^[^_]*_", "", cols), sub("_.*$", "", cols))
get_encode_info <- function(acc) c(epigenome = NA, target = unname(target[acc]))   # replaces the portal query

# ---- 1. feature matrix ----
if (file.exists(output_matrix)) {
    message("== reusing feature matrix ", output_matrix, " ==")
    fm <- readRDS(output_matrix)
} else {
    message("== building feature matrix ==")
    fm <- make_feature_matrix_from_bigwigs_and_bed(
        bigwig_dir = bigwig_dir, genomic_regions_bed = input_bed,
        exclude_files = exclude_files, num_cores = num_cores)
    saveRDS(fm, output_matrix)
}
message(sprintf("feature matrix: %d regions x %d features", nrow(fm), ncol(fm)))
if (!setequal(colnames(fm), cols))
    warning("columns differ from --columns: ", paste(setdiff(union(colnames(fm), cols), intersect(colnames(fm), cols)), collapse = ", "))

# ---- 2. neighbours of the methyldriver regions ----
md <- fread(md_path)
stopifnot(all(md$coordID %in% rownames(fm)))
x <- fm[md$coordID, , drop = FALSE]; rm(fm); invisible(gc())
rownames(x) <- md$mdID
message(sprintf("== finding most-similar regions: %d regions, %d neighbours, skipping < %d bp ==", nrow(x), n_neighbors, min_dist))
nb <- most_similar_n_regions_min_dist(x, md$chr, md$start, md$end, n = n_neighbors, n_pcs = n_pcs,
                                      min_dist = min_dist, num_cores = num_cores)
saveRDS(nb, output_neighbors)
message("done")
