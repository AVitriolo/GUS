#!/usr/bin/env Rscript
# =============================================================================
# batch_correct_fsva.R — LEAKAGE-SAFE RNA batch correction (frozen SVA)
#
# sva() on TRAIN builds the surrogate variables; fsva() projects TEST onto
# the frozen train SV basis. Corrected train + test written; QC PCA (fit on
# train, test projected) before vs after.
#
# mod = mod0 = ~1  -> protect nothing: the SVs absorb the dominant residual
# axis (CPGEA/MCRPC dataset). Because dataset is confounded with metastasis,
# removing it also removes metastasis biology. N<->T is the clean contrast.
#
# Parker, Corrada Bravo & Leek (2014), PeerJ 2:e561.
# =============================================================================
suppressPackageStartupMessages({
  library(sva); library(ggplot2)
})

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
counts_train_path <- args$counts_train
counts_test_path  <- args$counts_test
out_train         <- args$out_train
out_test          <- args$out_test
out_fit           <- args$out_fit
out_plot_cond     <- args$out_plot_cond
out_plot_ds       <- args$out_plot_ds
n_sv_arg          <- args$n_sv

train <- as.matrix(read.table(counts_train_path, header = TRUE,
                              row.names = 1, sep = "\t", check.names = FALSE))
test  <- as.matrix(read.table(counts_test_path,  header = TRUE,
                              row.names = 1, sep = "\t", check.names = FALSE))
stopifnot(identical(rownames(train), rownames(test)))
cat(sprintf("Train: %d genes x %d samples | Test: %d genes x %d samples\n",
            nrow(train), ncol(train), nrow(test), ncol(test)))

# ---- labels (for QC only; NOT protected in the model) ----
cond_of <- function(s) ifelse(grepl("^N", s), "Normal",
                       ifelse(grepl("^T", s), "Tumour",
                       ifelse(grepl("^DTB", s), "metastasis", NA)))
ds_of   <- function(s) ifelse(grepl("^DTB", s), "MCRPC", "CPGEA")

cond_train <- factor(cond_of(colnames(train)), levels = c("Normal","Tumour","metastasis"))
cond_test  <- factor(cond_of(colnames(test)),  levels = c("Normal","Tumour","metastasis"))
ds_train   <- factor(ds_of(colnames(train)),   levels = c("CPGEA","MCRPC"))
ds_test    <- factor(ds_of(colnames(test)),    levels = c("CPGEA","MCRPC"))

meta_train <- data.frame(
    condition = cond_train,
    dataset = ds_train,
    row.names = colnames(train)
)

# ---- 1. SVA on TRAIN only (protect nothing: mod = mod0 = ~1) ----
mod  <- model.matrix(~ 1, data = data.frame(cond_train))
mod0 <- model.matrix(~ 1, data = data.frame(cond_train))



cat(sprintf("Train NA: %d, NaN: %d, Inf: %d\n",
            sum(is.na(train)), sum(is.nan(train)), sum(is.infinite(train))))
cat(sprintf("Train range: [%.3f, %.3f]\n", min(train, na.rm=TRUE), max(train, na.rm=TRUE)))
# zero-variance genes break SVA's density estimation
rv <- apply(train, 1, var)
cat(sprintf("Zero-variance genes in train: %d\n", sum(rv == 0, na.rm=TRUE)))
cat(sprintf("Genes with NA variance: %d\n", sum(is.na(rv))))



n_sv <- if (!is.null(n_sv_arg)) as.integer(n_sv_arg) else num.sv(train, mod, method = "leek")
cat(sprintf("Estimating %d surrogate variable(s) on TRAIN (mod = ~1) ...\n", n_sv))
svobj <- sva(train, mod, mod0, n.sv = n_sv)
cat(sprintf("SVA returned %d SV(s).\n", svobj$n.sv))

# ---- 2. frozen SVA: adjust train DB + project/adjust TEST ----
cat("Running fsva (exact) ...\n")
fsvaobj <- fsva(dbdat = train, mod = mod, sv = svobj,
                newdat = test, method = "exact")
train_adj <- fsvaobj$db
test_adj  <- fsvaobj$new

write.table(train_adj, out_train, sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = TRUE)
write.table(test_adj,  out_test,  sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = TRUE)
saveRDS(list(svobj = svobj, newsv = fsvaobj$newsv, mod = mod, n_sv = n_sv),
        out_fit)

# ---- 3. QC PCA: fit on TRAIN, project TEST (before vs after) ----
.rowVars <- if (requireNamespace("matrixStats", quietly = TRUE))
  matrixStats::rowVars else function(m) apply(m, 1, var)

pca_pair <- function(mat_train, mat_test, stage, n_top = 2000) {
  vg  <- head(order(.rowVars(mat_train), decreasing = TRUE), min(n_top, nrow(mat_train)))
  Xtr <- t(mat_train[vg, , drop = FALSE]); ctr <- colMeans(Xtr)
  p   <- prcomp(Xtr, center = TRUE, scale. = FALSE)
  Xte <- scale(t(mat_test[vg, , drop = FALSE]), center = ctr, scale = FALSE) %*% p$rotation
  rbind(
    data.frame(PC1 = p$x[,1],  PC2 = p$x[,2],  condition = cond_train,
               dataset = ds_train, set = "train", stage = stage),
    data.frame(PC1 = Xte[,1],  PC2 = Xte[,2],  condition = cond_test,
               dataset = ds_test,  set = "test",  stage = stage))
}

df <- rbind(pca_pair(train, test, "before"),
            pca_pair(train_adj, test_adj, "after"))
df$stage <- factor(df$stage, levels = c("before","after"))
df$set   <- factor(df$set,   levels = c("train","test"))

cols_cond <- c(Normal="#2166ac", Tumour="#d6604d", metastasis="#1a9641")
cols_ds   <- c(CPGEA="#4575b4", MCRPC="#d73027")

g_cond <- ggplot(df, aes(PC1, PC2, color = condition)) +
  geom_point(size = 1.8, alpha = 0.8) + facet_grid(set ~ stage, scales = "free") +
  scale_color_manual(values = cols_cond) + theme_bw(base_size = 12) +
  labs(title = "RNA fsva — PCA by CONDITION",
       subtitle = "PCA fit on train; test projected")
g_ds <- ggplot(df, aes(PC1, PC2, color = dataset)) +
  geom_point(size = 1.8, alpha = 0.8) + facet_grid(set ~ stage, scales = "free") +
  scale_color_manual(values = cols_ds) + theme_bw(base_size = 12) +
  labs(title = "RNA fsva — PCA by DATASET (want: separation shrinks after)",
       subtitle = "PCA fit on train; test projected")

ggsave(out_plot_cond, g_cond, width = 10, height = 7)
ggsave(out_plot_ds,   g_ds,   width = 10, height = 7)

cat(sprintf("Saved:\n  %s\n  %s\n  %s\n  %s\n", out_train, out_test, out_plot_cond, out_plot_ds))