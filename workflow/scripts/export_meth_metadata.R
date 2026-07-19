#!/usr/bin/env Rscript
# ============================================================
# export_meth_metadata.R
# Export colData + rowRanges-chromosome as TSV, so the Python
# ComBat step (combat_correction_meth.py) can read the sample
# metadata (dataset/condition/split) and the per-CpG chromosome.
# ============================================================
suppressPackageStartupMessages({ library(GenomicRanges) })

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

cd <- readRDS(args$in_coldata)
cd$sample <- rownames(cd)
write.table(cd[, c("sample", "condition", "dataset", "split")],
            args$out_coldata, sep = "\t", quote = FALSE,
            row.names = TRUE, col.names = NA)

rr <- readRDS(args$in_rowranges)
write.table(data.frame(seqnames = as.character(seqnames(rr))),
            args$out_rowranges, sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("colData: %d samples -> %s\n", nrow(cd), args$out_coldata))
cat(sprintf("rowRanges: %d CpGs -> %s\n", length(rr), args$out_rowranges))