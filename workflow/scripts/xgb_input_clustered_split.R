options(scipen=999)

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

input_path_rse          <- args$input_path_rse
input_path_mapping      <- args$input_path_mapping
input_path_counts_train <- args$input_path_counts_train     # NEW
input_path_counts_test  <- args$input_path_counts_test      # NEW
train_samples_path      <- args$train_samples               # NEW
test_samples_path       <- args$test_samples                # NEW
output_path_train       <- args$output_path_train           # NEW
output_path_test        <- args$output_path_test            # NEW

TRAIN <- readLines(train_samples_path)
TEST  <- readLines(test_samples_path)

parts <- unlist(strsplit(basename(input_path_rse), split = "_"))
TxID <- parts[length(parts)]

# ---- load methylation matrix from rse ----
rse  <- HDF5Array::loadHDF5SummarizedExperiment(input_path_rse)
meth <- SummarizedExperiment::assay(rse, "M")
meth <- as.matrix(meth)
# meth: CpGs x samples

# ---- load cluster mapping (the FIT: derived from TRAIN correlations/expression) ----
mapping <- read.table(input_path_mapping, header = FALSE, sep = "\t")
colnames(mapping) <- c("CpGID", "cluster_id")
mapping$cluster_id <- paste0(TxID, "_", mapping$cluster_id)

# ---- APPLY the train-derived map: mean methylation per cluster, per sample ----
# All samples are averaged here; test samples never influenced WHICH CpGs are
# grouped together - they are simply averaged according to the train grouping.
common_cpgs <- intersect(mapping$CpGID, rownames(meth))
meth        <- meth[common_cpgs, , drop = FALSE]
mapping     <- mapping[mapping$CpGID %in% common_cpgs, ]

cluster_means <- sapply(split(mapping$CpGID, mapping$cluster_id), function(cpg_ids) {
    colMeans(meth[cpg_ids, , drop = FALSE], na.rm = TRUE)
})
# cluster_means: samples x clusters
cluster_means <- as.data.frame(cluster_means)

cat(sprintf("[%s] %d CpGs -> %d clusters, %d samples\n",
            TxID, length(common_cpgs), ncol(cluster_means), nrow(cluster_means)))

# ---- split into train / test, attaching the matching counts ----
build_set <- function(counts_path, samples, label, out_path) {

    counts <- read.table(counts_path, header = TRUE, sep = "\t",
                         row.names = 1, check.names = FALSE)
    if (!TxID %in% rownames(counts))
        stop(sprintf("[%s] not found in %s counts.", TxID, label))

    expr <- as.numeric(counts[TxID, ])
    names(expr) <- colnames(counts)

    s <- Reduce(intersect, list(rownames(cluster_means), names(expr), samples))
    if (!length(s))
        stop(sprintf("[%s] no %s samples shared between RSE and counts.", TxID, label))

    out <- cluster_means[s, , drop = FALSE]
    out[[TxID]] <- expr[s]

    cat(sprintf("[%s] %-5s: %d samples x %d metaCpGs\n",
                TxID, label, nrow(out), ncol(out) - 1))

    write.table(out, out_path, sep = "\t",
                row.names = TRUE, quote = FALSE, col.names = NA)
    out
}

xgb_train <- build_set(input_path_counts_train, TRAIN, "train", output_path_train)
xgb_test  <- build_set(input_path_counts_test,  TEST,  "test",  output_path_test)

# columns must match: the model trains on train's columns and predicts on test's
stopifnot(identical(colnames(xgb_train), colnames(xgb_test)))