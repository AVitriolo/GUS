# =============================================================================
# merge_benign_M_strands.R
# WCDT benign-adjacent WGBS (Bismark CpG reports, one row per CpG per strand)
# -> one row per CpG, positioned like the MCRPC tumour WGBS (the + strand C):
# minus-strand rows move to pos - 1 and counts are summed. Verified on
# NT-DTB-112-BL chr20: offset -1 maps 100% of CpGs onto MCRPC CpGs (0/+1: ~50%).
# Output per sample: chrom, pos, M, U (gzipped tsv), main chromosomes only.
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(parallel) })

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
in_dir  <- args$in_dir
out_dir <- args$out_dir
n_cores <- as.integer(if (is.null(args$n_cores)) 5 else args$n_cores)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

CHROMS <- paste0("chr", c(1:22, "X", "Y"))
files  <- list.files(in_dir, "\\.tsv\\.gz$", full.names = TRUE)

res <- mclapply(files, function(f) {
    x <- fread(f, header = FALSE, select = 1:5, col.names = c("chrom", "pos", "strand", "M", "U"))
    x <- x[chrom %in% CHROMS]
    x[strand == "-", pos := pos - 1L]
    y <- x[, .(M = sum(M), U = sum(U)), by = .(chrom, pos)]
    setorder(y, chrom, pos)
    out <- file.path(out_dir, sub("\\.tsv\\.gz$", ".merged.tsv.gz", basename(f)))
    fwrite(y, out, sep = "\t")
    sprintf("%s: %d CpGs, median cov %d", basename(f), nrow(y), as.integer(median(y$M + y$U)))
}, mc.cores = n_cores)
cat(unlist(res), sep = "\n")
