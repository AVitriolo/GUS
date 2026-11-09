options(scipen=999)

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

input_path_rse    <- args$input_path_rse
output_path_corr  <- args$output_path_corr

`%>%` <- magrittr::`%>%`

# ---- load RSE ----
rse  <- HDF5Array::loadHDF5SummarizedExperiment(input_path_rse)
meth <- SummarizedExperiment::assay(rse, "M")
meth <- as.matrix(meth)

CpGs <- SummarizedExperiment::rowRanges(rse)
coord_order <- GenomicRanges::mcols(GenomicRanges::sort(CpGs))$CpGID

# from here samples x CpGs
meth_t <- as.data.frame(t(meth))
colnames(meth_t) <- rownames(meth)

# order columns by genomic coordinate
coord_order <- coord_order[coord_order %in% colnames(meth_t)]
meth_t <- meth_t[, coord_order, drop = FALSE]

correlations_table <- cor(meth_t, method = "spearman", use = "pairwise.complete.obs")

write.table(correlations_table,
            file      = output_path_corr,
            col.names = TRUE,
            row.names = TRUE,
            quote     = FALSE,
            sep       = "\t")