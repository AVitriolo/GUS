options(scipen=999)

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

input_path_rse     <- args$input_path_rse
input_path_mapping <- args$input_path_mapping
input_path_counts  <- args$input_path_counts
output_path        <- args$output_path

parts <- unlist(strsplit(basename(input_path_rse), split = "_"))
TxID <- parts[length(parts)]

# ---- load methylation matrix from rse ----
rse  <- HDF5Array::loadHDF5SummarizedExperiment(input_path_rse)
meth <- SummarizedExperiment::assay(rse, "M")
meth <- as.matrix(meth)
# meth: CpGs x samples

# ---- load cluster mapping ----
mapping <- read.table(input_path_mapping, header = FALSE, sep = "\t")
colnames(mapping) <- c("CpGID", "cluster_id")
mapping$cluster_id <- paste0(TxID, "_", mapping$cluster_id)

# ---- compute mean methylation per cluster, per sample ----
common_cpgs <- intersect(mapping$CpGID, rownames(meth))
meth <- meth[common_cpgs, , drop = FALSE]
mapping <- mapping[mapping$CpGID %in% common_cpgs, ]

cluster_means <- sapply(split(mapping$CpGID, mapping$cluster_id), function(cpg_ids) {
    colMeans(meth[cpg_ids, , drop = FALSE], na.rm = TRUE)
})
# cluster_means: samples x clusters

cluster_means <- as.data.frame(cluster_means)

# ---- load expression for this TxID ----
counts <- read.table(input_path_counts, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
expr   <- as.numeric(counts[TxID, ])
names(expr) <- colnames(counts)

# ---- match samples ----
common_samples <- intersect(rownames(cluster_means), names(expr))
cluster_means <- cluster_means[common_samples, , drop = FALSE]
cluster_means[[TxID]] <- expr[common_samples]

write.table(cluster_means, output_path, sep = "\t", row.names = TRUE, quote = FALSE, col.names = NA)
