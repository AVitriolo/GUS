#!/usr/bin/env Rscript
# =============================================================================
# assemble_featmatrix_w20.R
# Regions x bigWig matrix from the per-file averages of
# bigwig_stream_average.sh (name, mean), rows matched by region name.
# Columns named <accession>_<target> as in make_feature_matrix_from_bigwigs_and_bed.R;
# the targets come from --columns (the column names of the HPC-built matrix),
# since the ENCODE portal does not answer from this machine.
# =============================================================================
suppressPackageStartupMessages(library(data.table))

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
bed      <- fread(args$bed, header = FALSE, select = 4, col.names = "name")$name
cols     <- readLines(args$columns)                     # ENCFF..._<target>
part_dir <- args$part_dir
output   <- args$output

acc   <- sub("_.*$", "", cols)
files <- file.path(part_dir, paste0(acc, ".tsv.gz"))
if (any(!file.exists(files))) stop("missing averages for: ", paste(acc[!file.exists(files)], collapse = ", "))

fm <- matrix(NA_real_, nrow = length(bed), ncol = length(cols), dimnames = list(bed, cols))
for (j in seq_along(files)) {
    d <- fread(files[j], header = FALSE, col.names = c("name", "mean"))
    i <- match(bed, d$name)
    if (anyNA(i)) stop(sprintf("%s: %d regions missing", acc[j], sum(is.na(i))))
    fm[, j] <- d$mean[i]
    if (j %% 100 == 0) cat(sprintf("%d / %d\n", j, length(files)))
}
cat(sprintf("feature matrix: %d regions x %d features\n", nrow(fm), ncol(fm)))
saveRDS(fm, output)
