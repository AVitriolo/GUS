# =============================================================================
# purity_sensitivity_vs_performance.R
# Two checks of the purity correction that use the real XGBoost results only
# (per-transcript performance + selected CpGs), no new models.
#  1) Is the bulk -> purified change in R2 larger for purity-sensitive
#     transcripts (expression correlated with purity in tumours)?
#  2) Do the models select purity-proxy CpGs, and less so after correction?
# Purity association = pooled within-group Pearson r (T and M centred
# separately), so the T-vs-M purity/batch difference cannot drive it.
# =============================================================================
suppressPackageStartupMessages({
    library(HDF5Array); library(GenomicRanges); library(data.table); library(parallel)
})
args   <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
ncore  <- as.integer(if (is.null(args$ncore)) 4 else args$ncore)
outdir <- if (is.null(args$outdir)) "results/purity_checks/purity_sensitivity" else args$outdir
dir.create(file.path(outdir, "cpg_r"), recursive = TRUE, showWarnings = FALSE)
PUR     <- "/home/ieo7293/Prostate_GAS_purified"
BED_DIR <- "resources/CpG_data/top_K_2Mb_full"
PERF_B  <- "results/performance_2mb_full_kcapped_cv3r5"
SEL_B   <- "results/sel_feats_2mb_full_kcapped_cv3r5"
PERF_P  <- file.path(PUR, "results/performance")
SEL_P   <- file.path(PUR, "results/sel_feats")

# ---- purity (PAMES), tumours only -------------------------------------------
pn <- read.delim("resources/purity/pames_purity_NT.tsv"); pm <- read.delim("resources/purity/pames_purity_M.tsv")
pn <- pn[pn$Condition == "Tumour", ]
pur <- c(setNames(pn$Purity, pn$Sample), setNames(pm[[2]], pm[[1]]))
grp_of <- function(s) ifelse(grepl("^DTB", s), "M", "T")

# pooled within-group r of every row of x (features x samples) with p
wg_cor <- function(x, p, g) {
    for (k in unique(g)) { j <- g == k; x[, j] <- x[, j] - rowMeans(x[, j, drop = FALSE]); p[j] <- p[j] - mean(p[j]) }
    as.numeric((x %*% p) / sqrt(rowSums(x^2) * sum(p^2)))
}

# ---- 1) CpG x purity, both stores (per-chr checkpoints) ----------------------
gr <- readRDS("resources/integrated_wgbs_hdf5/integrated_wgbs_combat_rowRanges.rds")
stores <- list(
    bulk = list(h5 = "resources/integrated_wgbs_hdf5/integrated_wgbs_combat.h5",
                s  = rownames(readRDS("resources/integrated_wgbs_hdf5/integrated_wgbs_combat_colData.rds"))),
    pur  = list(h5 = file.path(PUR, "resources/integrated_wgbs_hdf5/integrated_wgbs_combat.h5"),
                s  = rownames(readRDS(file.path(PUR, "resources/integrated_wgbs_hdf5/integrated_wgbs_combat_colData.rds")))))
chr_rows <- split(seq_along(gr), as.character(seqnames(gr)))
# same tumours in both stores, so the two r's are comparable
tum <- Reduce(intersect, list(names(pur), stores$bulk$s, stores$pur$s))
cat(sprintf("tumours with PAMES purity in both WGBS stores: %d (T %d, M %d)\n",
            length(tum), sum(grp_of(tum) == "T"), sum(grp_of(tum) == "M")))

invisible(mclapply(names(chr_rows), function(ch) {
    f <- file.path(outdir, "cpg_r", paste0(ch, ".rds"))
    if (file.exists(f)) return(NULL)
    rows <- chr_rows[[ch]]
    res <- lapply(names(stores), function(st) {
        M <- HDF5Array(stores[[st]]$h5, "M"); j <- match(tum, stores[[st]]$s)
        blocks <- split(seq_along(rows), ceiling(seq_along(rows) / 5e5))
        unlist(lapply(blocks, function(b) {
            x <- as.matrix(M[rows[b], j, drop = FALSE])
            r <- wg_cor(x, pur[tum], grp_of(tum)); r[!is.finite(r)] <- NA; r
        }), use.names = FALSE)
    })
    out <- data.table(pos = start(gr)[rows], r_bulk = res[[1]], r_pur = res[[2]])
    saveRDS(out, paste0(f, ".tmp")); file.rename(paste0(f, ".tmp"), f)
    cat(sprintf("[%s] %s: %d CpGs\n", format(Sys.time(), "%H:%M"), ch, nrow(out)))
    NULL
}, mc.cores = ncore, mc.preschedule = FALSE))
cpg_r <- rbindlist(lapply(names(chr_rows), function(ch) {
    d <- readRDS(file.path(outdir, "cpg_r", paste0(ch, ".rds"))); d[, id := paste0(ch, "_", pos)]; d
}))
setkey(cpg_r, id)
fwrite(data.table(r_bulk = quantile(abs(cpg_r$r_bulk), c(.5, .9, .99), na.rm = TRUE),
                  r_pur  = quantile(abs(cpg_r$r_pur),  c(.5, .9, .99), na.rm = TRUE), q = c(.5, .9, .99)),
       file.path(outdir, "cpg_abs_r_quantiles.tsv"), sep = "\t")

# ---- 2) transcript expression x purity ---------------------------------------
rd <- function(f) as.matrix(read.delim(f, row.names = 1, check.names = FALSE))
expr_r <- function(dir) {
    x <- cbind(rd(file.path(dir, "resources/RNA/kallisto_counts_NTM_train_combat_25_10.tsv")),
               rd(file.path(dir, "resources/RNA/kallisto_counts_NTM_test_combat_25_10.tsv")))
    s <- intersect(colnames(x), names(pur))
    list(r = setNames(wg_cor(x[, s], pur[s], grp_of(s)), rownames(x)), n = length(s))
}
eB <- expr_r("."); eP <- expr_r(PUR)
cat(sprintf("tumours with purity in expression: bulk %d, purified %d\n", eB$n, eP$n))

# ---- 3) per transcript: performance + purity-proxy score of selected CpGs ----
kv <- function(f) { d <- fread(f, header = FALSE); setNames(d$V2, d$V1) }
tx_of <- function(f) sub(".*_(ENST[0-9]+)$", "\\1", f)
fb <- list.files(PERF_B); fp <- list.files(PERF_P, pattern = "_1000000_")
txs <- intersect(tx_of(fb), tx_of(fp))
fb <- setNames(fb, tx_of(fb)); fp <- setNames(fp, tx_of(fp))
sb <- list.files(SEL_B); sb <- setNames(sb, tx_of(sb))
sp <- list.files(SEL_P, pattern = "_1000000_"); sp <- setNames(sp, tx_of(sp))
beds <- list.files(BED_DIR); beds <- setNames(beds, sub(".*_(ENST[0-9]+)\\.bed$", "\\1", beds))
cat(sprintf("transcripts with bulk + purified performance: %d\n", length(txs)))

mabs <- function(ids, col) { v <- cpg_r[J(ids), get(col)]; mean(abs(v), na.rm = TRUE) }
sel_ids <- function(dir, f) if (is.na(f)) character(0) else sub("^ENST[0-9]+_", "", readLines(file.path(dir, f)))

tab <- rbindlist(mclapply(txs, function(tx) {
    a <- kv(file.path(PERF_B, fb[[tx]])); b <- kv(file.path(PERF_P, fp[[tx]]))
    ib <- sel_ids(SEL_B, sb[tx]); ip <- sel_ids(SEL_P, sp[tx])
    iw <- if (is.na(beds[tx])) character(0) else fread(file.path(BED_DIR, beds[[tx]]), header = FALSE)$V2
    data.table(tx = tx,
               cv_bulk = a[["cv_r2"]], cv_pur = b[["cv_r2"]],
               test_bulk = a[["test_r2"]], test_pur = b[["test_r2"]],
               r_expr_bulk = unname(eB$r[tx]), r_expr_pur = unname(eP$r[tx]),
               k_bulk = length(ib), k_pur = length(ip), n_window = length(iw),
               # purity-proxy score = mean |r(CpG, purity)|, measured in BULK methylation
               win_absr = mabs(iw, "r_bulk"), selB_absr = mabs(ib, "r_bulk"), selP_absr = mabs(ip, "r_bulk"),
               # the same, measured in the data each model was trained on
               selP_absr_own = mabs(ip, "r_pur"), win_absr_pur = mabs(iw, "r_pur"))
}, mc.cores = ncore))
fwrite(tab, file.path(outdir, "purity_sensitivity_per_tx.tsv.gz"), sep = "\t")
cat(sprintf("written %d transcripts -> %s\n", nrow(tab), outdir))
