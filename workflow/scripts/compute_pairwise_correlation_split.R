options(scipen=999)

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

input_path_rse     <- args$input_path_rse
train_samples_path <- args$train_samples         
output_path_corr   <- args$output_path_corr
TxID               <- args$tx_id
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

# ---- degenerate case: a single CpG cannot be correlated pairwise ----
# cor() on a 1-column matrix errors ("non-conformable arrays"). A single CpG
# trivially forms one cluster, so a 1x1 correlation of 1 is the correct output.

if (ncol(meth_t) < 2) {
  cat(sprintf("[%s] only %d CpG(s) after train subset - writing trivial 1x1 correlation\n",
              TxID, ncol(meth_t)))
  corr <- matrix(1, nrow = ncol(meth_t), ncol = ncol(meth_t),
                 dimnames = list(colnames(meth_t), colnames(meth_t)))
  write.table(corr,
              file      = output_path_corr,
              col.names = TRUE,
              row.names = TRUE,
              quote     = FALSE,
              sep       = "\t")
  cat(sprintf("[%s] written %d x %d correlation matrix\n",
              TxID, nrow(corr), ncol(corr)))
  quit(save = "no", status = 0)
}

# order columns by genomic coordinate (only reached when >= 2 CpGs)
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