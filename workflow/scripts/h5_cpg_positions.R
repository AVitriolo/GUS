# chr / pos of every CpG in the integrated h5, in h5 column order (input of the stability-selection scripts)
# Usage: Rscript h5_cpg_positions.R <rowRanges.rds> <out.tsv>
suppressPackageStartupMessages(library(GenomicRanges))
args <- commandArgs(trailingOnly = TRUE)
rr <- readRDS(args[1])
out <- data.frame(chr = as.character(seqnames(rr)), pos = start(rr))
write.table(out, args[2], sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("%d CpGs -> %s\n", nrow(out), args[2]))
