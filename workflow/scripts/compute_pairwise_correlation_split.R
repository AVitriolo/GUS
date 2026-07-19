options(scipen=999)

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

input_path_rse     <- args$input_path_rse
train_samples_path <- args$train_samples         
output_path_corr   <- args$output_path_corr

`%>%` <- magrittr::`%>%`

TRAIN <- readLines(train_samples_path)

# ---- load RSE ----
rse  <- HDF5Array::loadHDF5SummarizedExperiment(input_path_rse)
meth <- SummarizedExperiment::assay(rse, "M")
meth <- as.matrix(meth)                            # CpGs x samples

CpGs <- SummarizedExperiment::rowRanges(rse)
coord_order <- GenomicRanges::mcols(GenomicRanges::sort(CpGs))$CpGID

# from here samples x CpGs
meth_t <- as.data.frame(t(meth))
colnames(meth_t) <- rownames(meth)

# ---- LEAKAGE-SAFE: restrict to TRAIN samples --------------------------------
n_before <- nrow(meth_t)
meth_t   <- meth_t[rownames(meth_t) %in% TRAIN, , drop = FALSE]
cat(sprintf("Samples: %d -> %d train\n", n_before, nrow(meth_t)))

if (nrow(meth_t) < 3)
  stop(sprintf("Only %d train samples - cannot compute correlations.", nrow(meth_t)))

# order columns by genomic coordinate
coord_order <- coord_order[coord_order %in% colnames(meth_t)]
meth_t <- meth_t[, coord_order, drop = FALSE]

correlations_table <- cor(meth_t, method = "spearman", use = "pairwise.complete.obs")

# CpGs with zero variance across TRAIN samples yield NA
n_na <- sum(is.na(correlations_table))
if (n_na > 0) {
  cat(sprintf("WARNING: %d NA correlations (zero-variance CpGs in train); set to 0\n", n_na))
  correlations_table[is.na(correlations_table)] <- 0
  diag(correlations_table) <- 1
}

write.table(correlations_table,
            file      = output_path_corr,
            col.names = TRUE,
            row.names = TRUE,
            quote     = FALSE,
            sep       = "\t")

cat(sprintf("Written %d x %d correlation matrix (from %d train samples)\n",
            nrow(correlations_table), ncol(correlations_table), nrow(meth_t)))