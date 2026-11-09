# Stability selection with the models' own selector, step 2: per-gene stability and per-CpG frequency.
# Input: subsample_selections_shard*.tsv / candidates_shard*.tsv from stability_selection_kbest.py
# (B = 2 x n_pairs complementary half-samples of the train samples, SelectKBest f_regression, gene's k).
#
# Per gene:
#   phi          Nogueira et al. (2018) stability index of the B selections (stabm::stabilityNogueira),
#                ~0 = no more stable than random selection of k CpGs, 1 = identical selections
#   cutoff_pfer1 frequency cutoff giving PFER <= 1 (Shah & Samworth complementary-pairs bound,
#                unimodality assumption; stabs::stabsel_parameters), NA if not attainable for this q/p
#   n_stable     CpGs with frequency >= cutoff_pfer1
#   pct_saved_stable  % of the model's saved selected CpGs with frequency >= cutoff_pfer1
# Per CpG (cpg_freq.tsv): selection frequency over the B subsamples, saved_selected flag, stable flag.
#
# Usage: Rscript stability_selection_summary.R <in_dir> [n_pairs=100] [pfer=1]

suppressPackageStartupMessages({ library(stabm); library(stabs); library(data.table) })

args    <- commandArgs(trailingOnly = TRUE)
in_dir  <- args[1]
n_pairs <- if (length(args) >= 2) as.integer(args[2]) else 100L
pfer    <- if (length(args) >= 3) as.numeric(args[3]) else 1

sel  <- rbindlist(lapply(list.files(in_dir, "^subsample_selections_shard.*\\.tsv$", full.names = TRUE),
                         fread, sep = "\t", header = FALSE, col.names = c("TxID", "b", "idx"), colClasses = "character"))
cand <- rbindlist(lapply(list.files(in_dir, "^candidates_shard.*\\.tsv$", full.names = TRUE),
                         fread, sep = "\t", header = FALSE,
                         col.names = c("TxID", "n_cand", "k", "idx", "CpGID", "saved_selected")))
B <- 2L * n_pairs

cutoff_for <- function(p, q) tryCatch(
    stabsel_parameters(p = p, q = q, PFER = pfer, B = n_pairs,
                       sampling.type = "SS", assumption = "unimodal")$cutoff,
    error = function(e) NA_real_)

info <- unique(cand[, .(TxID, n_cand, k)])
sel  <- merge(sel, info, by = "TxID")
genes <- sel[, .(n_sub = .N, phi = stabilityNogueira(lapply(strsplit(idx, ","), function(s) as.integer(s) + 1L),
                                                     p = n_cand[1])),
             by = .(TxID, k, n_cand)]
stopifnot(all(genes$n_sub == B))

freq <- sel[, .(idx = as.integer(unlist(strsplit(idx, ",")))), by = TxID][, .(freq = .N / B), by = .(TxID, idx)]
cand <- merge(cand, freq, by = c("TxID", "idx"), all.x = TRUE)
cand[is.na(freq), freq := 0]

genes[, cutoff_pfer1 := mapply(cutoff_for, n_cand, k)]
cand <- merge(cand, genes[, .(TxID, cutoff_pfer1)], by = "TxID")
cand[, stable := !is.na(cutoff_pfer1) & freq >= cutoff_pfer1]

genes <- merge(genes, cand[, .(n_stable = sum(stable),
                               pct_saved_stable = 100 * mean(stable[saved_selected == 1]),
                               mean_freq_saved = mean(freq[saved_selected == 1])), by = TxID], by = "TxID")

fwrite(genes, file.path(in_dir, "stability_genes.tsv"), sep = "\t")
fwrite(cand[order(TxID, -freq)], file.path(in_dir, "cpg_freq.tsv"), sep = "\t")
cat(sprintf("%d genes; phi median %.3f; cutoff NA (PFER<=%g not attainable) for %d genes\n",
            nrow(genes), median(genes$phi), pfer, sum(is.na(genes$cutoff_pfer1))))
print(summary(genes[, .(phi, n_stable, pct_saved_stable, mean_freq_saved)]))
