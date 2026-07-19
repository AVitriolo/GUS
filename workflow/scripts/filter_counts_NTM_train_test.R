source("workflow/scripts/helpers/filterGenes.R")
source("workflow/scripts/helpers/do_norm_fit_apply.R")
source("workflow/scripts/helpers/load_and_clean.R")

options(scipen = 999)
args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

input_expr_values_path_NT <- args$input_expr_values_path_NT
input_expr_values_path_M  <- args$input_expr_values_path_M
train_samples_path        <- args$train_samples_path
test_samples_path         <- args$test_samples_path
minCount_expr             <- as.numeric(args$minCount_expr)
minSamples_expr           <- as.numeric(args$minSamples_expr)
output_path_counts_train  <- args$output_path_counts_train
output_path_counts_test   <- args$output_path_counts_test
output_path_TxIDs         <- args$output_path_TxIDs
output_path_tmm_fit       <- args$output_path_tmm_fit

NT <- load_and_clean(input_expr_values_path_NT)
M  <- load_and_clean(input_expr_values_path_M)
raw <- cbind(NT, M[rownames(NT), , drop = FALSE])   # align on TxIDs

train_samples <- readLines(train_samples_path)
test_samples  <- readLines(test_samples_path)
stopifnot(all(train_samples %in% colnames(raw)), all(test_samples %in% colnames(raw)))

raw_train <- raw[, train_samples, drop = FALSE]
raw_test  <- raw[, test_samples,  drop = FALSE]
cat(sprintf("Train: %d samples | Test: %d samples\n", ncol(raw_train), ncol(raw_test)))

# ===== FIT ON TRAIN =====
# per-condition filtering (TRAIN ONLY)
N_tr <- raw_train[, grepl("^N",   colnames(raw_train)), drop = FALSE]
T_tr <- raw_train[, grepl("^T",   colnames(raw_train)), drop = FALSE]
M_tr <- raw_train[, grepl("^DTB", colnames(raw_train)), drop = FALSE]
cat(sprintf("Train samples - N: %d, T: %d, M: %d\n", ncol(N_tr), ncol(T_tr), ncol(M_tr)))

N_ff <- filterGenes(N_tr, minCount = minCount_expr, minSamples = minSamples_expr)
T_ff <- filterGenes(T_tr, minCount = minCount_expr, minSamples = minSamples_expr)
M_ff <- filterGenes(M_tr, minCount = minCount_expr, minSamples = minSamples_expr)

# intersect TxIDs across the 3 conditions (TRAIN-DERIVED GENE LIST)
common_txids <- Reduce(intersect, list(rownames(N_ff), rownames(T_ff), rownames(M_ff)))
cat(sprintf("TxIDs passing (train) - N: %d, T: %d, M: %d, common: %d\n",
            nrow(N_ff), nrow(T_ff), nrow(M_ff), length(common_txids)))

# normalize TRAIN (all conditions together) on the common gene list
train_counts <- raw_train[common_txids, , drop = FALSE]
fit <- donorm_fit(train_counts, method = "TMM")
cat(sprintf("TMM fit on train (reference sample: %s)\n", fit$ref_name))

# ===== APPLY TO TEST =====
test_counts <- raw_test[common_txids, , drop = FALSE]   # SAME gene list
test_norm   <- donorm_apply(test_counts, fit)           # SAME reference + constant

# ---- log-transform both ----
log_train <- log2(as.matrix(fit$normalized) + 1)
log_test  <- log2(as.matrix(test_norm)      + 1)

write.table(log_train, file = output_path_counts_train,
            col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")
write.table(log_test,  file = output_path_counts_test,
            col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")
writeLines(common_txids, con = output_path_TxIDs, sep = "\n")

saveRDS(fit, output_path_tmm_fit)  # store params for the record
cat("Done.\n")
