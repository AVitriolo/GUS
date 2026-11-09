# Builds the CpG list for the test-set replication analysis (test_replication_selected_cpgs.py):
# for each modelled gene, its SELECTED CpGs (sel_feats) plus an equal-size RANDOM baseline of
# non-selected CpGs drawn from the same gene's 2Mb candidate set, each mapped to its column in
# the integrated HDF5 (M dataset is samples x CpGs, columns in rowRanges order).
#
# Usage: Rscript test_replication_index.R <sel_feats_dir> <candidate_bed_dir> <rowRanges.rds> <out.tsv> [max_random=50]

suppressPackageStartupMessages(library(GenomicRanges))

args <- commandArgs(trailingOnly = TRUE)
sel_dir    <- args[1]
cand_dir   <- args[2]
rowranges  <- args[3]
out_path   <- args[4]
max_random <- if (length(args) >= 5) as.integer(args[5]) else 50L

set.seed(1)
files <- list.files(sel_dir, full.names = TRUE)
out <- lapply(files, function(f) {
    TxID <- sub(".*_", "", basename(f))
    sel  <- readLines(f)
    sel  <- sub(paste0("^", TxID, "_"), "", sel)
    if (!length(sel)) return(NULL)

    cand_file <- file.path(cand_dir, sub("_3_ENST", "_ENST", basename(f)))
    cand_file <- paste0(cand_file, ".bed")
    rnd <- character(0)
    if (file.exists(cand_file)) {
        cand <- setdiff(read.table(cand_file, sep = "\t", colClasses = "character")[[2]], sel)
        rnd  <- sample(cand, min(length(cand), length(sel), max_random))
    }
    data.frame(TxID = TxID,
               CpGID = c(sel, rnd),
               set = rep(c("selected", "random"), c(length(sel), length(rnd))))
})
out <- do.call(rbind, out)

rr <- readRDS(rowranges)
out$h5_col <- match(out$CpGID, paste0(as.character(seqnames(rr)), "_", start(rr))) - 1L   # 0-based for h5py
rm(rr); invisible(gc())

n_miss <- sum(is.na(out$h5_col))
if (n_miss) cat(sprintf("WARNING: %d CpGs not found in rowRanges, dropped\n", n_miss))
out <- out[!is.na(out$h5_col), ]

cat(sprintf("%d genes, %d selected, %d random CpG-gene pairs\n",
            length(unique(out$TxID)), sum(out$set == "selected"), sum(out$set == "random")))
write.table(out, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
