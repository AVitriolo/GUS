# =============================================================================
# test_bmind_vs_demixt_NT.R
# Can reference-free bMIND (MIND package, Wang et al. Nat Commun 2021) replace
# DeMixT for purifying Tumour RNA? Tested on CPGEA train Tumours, where the
# DeMixT result with a real normal reference (PAMES purity) exists.
#
#   bMIND input : log2(CPM + 1) of the Tumours (CPM over ALL transcripts of
#                 each sample), fractions Tumour = PAMES, Other = 1 - PAMES,
#                 profile = covariance = NULL (no reference: prior from data)
#   subset      : marker transcripts + a random sample of pipeline transcripts
#
# Compared per transcript with DeMixT (demixt_purified_tumour_train_givenpi):
#   - correlation of per-sample purified values across Tumours
#   - agreement of the purification shift (purified - bulk, log2 scale)
#   - how much of the bulk tumour-to-tumour spread each method keeps
# =============================================================================
suppressPackageStartupMessages({ library(MIND); library(data.table) })

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
counts_path <- args$counts
tx_path     <- args$transcripts
cd_path     <- args$coldata
pames_path  <- args$pames
conv_path   <- args$conversion
demixt_dir  <- args$demixt_dir
n_random    <- as.integer(if (is.null(args$n_random)) 1500 else args$n_random)
ncore       <- as.integer(if (is.null(args$ncore)) 16 else args$ncore)
outdir      <- args$outdir
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

MARKER_GENES <- c("KRT19", "EPCAM", "CDH1", "PTPRC", "CD4", "CD163", "VWF", "ENG", "CDH5", "ACTA2", "DCN", "FAP",
                  "MYH11", "COL1A1", "COL3A1", "AMACR", "PCA3", "HPN", "FOLH1", "TDRD1", "ERG")

cd <- as.data.frame(readRDS(cd_path))
tumours <- rownames(cd)[cd$condition == "Tumour" & cd$split == "train"]
pm <- read.delim(pames_path, row.names = 1)
pi <- setNames(pm[[grep("^purity$", colnames(pm), ignore.case = TRUE)]], rownames(pm))[tumours]

raw <- fread(counts_path)
m   <- as.matrix(raw[, ..tumours]); rownames(m) <- raw$transcript_id
lcpm <- log2(t(t(m) / colSums(m)) * 1e6 + 1)          # CPM over the whole transcriptome

# transcripts: markers in the DeMixT output + random pipeline transcripts DeMixT kept
dx_P <- fread(file.path(demixt_dir, "demixt_purified_tumour_train_givenpi.tsv.gz"))
dx_B <- fread(file.path(demixt_dir, "bulk_scaled_plus1_tumour_train_givenpi.tsv.gz"))
dx_B <- dx_B[match(dx_P$transcript_id, dx_B$transcript_id)]
pipe <- fread(tx_path, select = 1)[[1]]
cv   <- suppressWarnings(fread(conv_path))[, .(transcript_id, gene_name)]
cand <- intersect(intersect(dx_P$transcript_id, pipe), rownames(lcpm))
cand <- cand[apply(lcpm[cand, , drop = FALSE], 1, sd) > 0]
mk   <- cand[cv$gene_name[match(cand, cv$transcript_id)] %in% MARKER_GENES]
set.seed(1)
tx   <- unique(c(mk, sample(setdiff(cand, mk), n_random)))
X    <- lcpm[tx, ]
cat(sprintf("%d train Tumours | %d transcripts (%d marker transcripts) | max log2CPM %.1f\n",
            ncol(X), nrow(X), length(mk), max(X)))

frac <- data.frame(Tumour = pi, Other = 1 - pi, row.names = tumours)
t0 <- Sys.time()
fit <- bMIND(X, frac = frac, ncore = ncore)
mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
A <- fit$A[, "Tumour", ]                               # transcripts x samples (sorted)
A <- A[, tumours]
cat(sprintf("bMIND: %d of %d transcripts returned in %.1f min (%.2f s per transcript on %d cores)\n",
            nrow(A), nrow(X), mins, 60 * mins / nrow(X), ncore))
saveRDS(fit, file.path(outdir, "bmind_fit_train_subset.rds"))

# ---- compare with DeMixT, per transcript ----
g   <- rownames(A)
dxP <- as.matrix(dx_P[match(g, transcript_id), ..tumours])
dxB <- as.matrix(dx_B[match(g, transcript_id), ..tumours])
rownames(dxP) <- rownames(dxB) <- g
shift_bm <- A - X[g, ]                                 # log2 scale
shift_dx <- log2(dxP + 1) - log2(dxB + 1)

per_tx <- data.table(
    transcript_id = g, gene = cv$gene_name[match(g, cv$transcript_id)],
    cor_purified = sapply(g, function(i) suppressWarnings(cor(A[i, ], log2(dxP[i, ] + 1), method = "spearman"))),
    median_shift_bmind = apply(shift_bm, 1, median), median_shift_demixt = apply(shift_dx, 1, median),
    spread_kept_bmind  = apply(A, 1, sd) / apply(X[g, ], 1, sd),
    spread_kept_demixt = apply(log2(dxP + 1), 1, sd) / apply(log2(dxB + 1), 1, sd))
fwrite(per_tx, file.path(outdir, "bmind_vs_demixt_per_transcript.tsv"), sep = "\t")

cat(sprintf("per-transcript Spearman (bMIND vs DeMixT purified, across Tumours): median %.2f, IQR %.2f-%.2f\n",
            median(per_tx$cor_purified, na.rm = TRUE), quantile(per_tx$cor_purified, .25, na.rm = TRUE),
            quantile(per_tx$cor_purified, .75, na.rm = TRUE)))
cat(sprintf("purification shift agreement across transcripts: Pearson %.2f | same sign %.0f%%\n",
            cor(per_tx$median_shift_bmind, per_tx$median_shift_demixt),
            100 * mean(sign(per_tx$median_shift_bmind) == sign(per_tx$median_shift_demixt))))
cat(sprintf("tumour-to-tumour spread kept (SD purified / SD bulk): bMIND median %.2f | DeMixT median %.2f\n",
            median(per_tx$spread_kept_bmind, na.rm = TRUE), median(per_tx$spread_kept_demixt, na.rm = TRUE)))

mk_tab <- per_tx[gene %in% MARKER_GENES, .(n_tx = .N,
    shift_bmind = round(median(median_shift_bmind), 2), shift_demixt = round(median(median_shift_demixt), 2),
    cor_purified = round(median(cor_purified, na.rm = TRUE), 2)), by = gene][order(match(gene, MARKER_GENES))]
print(mk_tab)
fwrite(mk_tab, file.path(outdir, "bmind_vs_demixt_markers.tsv"), sep = "\t")
