# =============================================================================
# residual_vs_purity_lambda.R
# Mild purity correction by blending: X_l = l*corrected + (1-l)*bulk for BOTH
# methylation (CAMDAC store vs bulk store) and expression (bMIND vs bulk),
# l in a grid. Per l: per-sample test residuals and per-transcript test R2.
# Same CpGs / samples / split for every l; only the blend weight changes.
# =============================================================================
suppressPackageStartupMessages({
    library(HDF5Array); library(GenomicRanges); library(data.table); library(parallel)
})
args    <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
K       <- as.integer(if (is.null(args$k)) 30 else args$k)
ncore   <- as.integer(if (is.null(args$ncore)) 16 else args$ncore)
LAMBDA  <- as.numeric(strsplit(if (is.null(args$lambda)) "0,0.25,0.5,0.75,1" else args$lambda, ",")[[1]])
bed_dir <- args$bed_dir; outdir <- args$outdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
PUR <- "/home/ieo7293/Prostate_GAS_purified"

gr <- readRDS("resources/integrated_wgbs_hdf5/integrated_wgbs_combat_rowRanges.rds")
sA <- rownames(readRDS("resources/integrated_wgbs_hdf5/integrated_wgbs_combat_colData.rds"))
sB <- rownames(readRDS(file.path(PUR, "resources/integrated_wgbs_hdf5/integrated_wgbs_combat_colData.rds")))
MB <- HDF5Array::HDF5Array("resources/integrated_wgbs_hdf5/integrated_wgbs_combat.h5", "M")
MA <- HDF5Array::HDF5Array(file.path(PUR, "resources/integrated_wgbs_hdf5/integrated_wgbs_combat.h5"), "M")

rd <- function(f) as.matrix(read.delim(f, row.names = 1, check.names = FALSE))
trB <- rd("resources/RNA/kallisto_counts_NTM_train_combat_25_10.tsv")
teB <- rd("resources/RNA/kallisto_counts_NTM_test_combat_25_10.tsv")
trA <- rd(file.path(PUR, "resources/RNA/kallisto_counts_NTM_train_combat_25_10.tsv"))
teA <- rd(file.path(PUR, "resources/RNA/kallisto_counts_NTM_test_combat_25_10.tsv"))

s_tr <- Reduce(intersect, list(colnames(trA), colnames(trB), sA, sB))
s_te <- Reduce(intersect, list(colnames(teA), colnames(teB), sA, sB))
keep <- c(s_tr, s_te); iA <- match(keep, sA); iB <- match(keep, sB)
jtr <- seq_along(s_tr); jte <- length(s_tr) + seq_along(s_te)
txs <- Reduce(intersect, list(rownames(trA), rownames(trB), rownames(teA), rownames(teB)))
cat(sprintf("samples: %d train, %d test | expression transcripts in both: %d | lambda: %s\n",
            length(s_tr), length(s_te), length(txs), paste(LAMBDA, collapse = ",")))

beds  <- list.files(bed_dir, full.names = TRUE)
tx_of <- sub(".*_(ENST[0-9]+)\\.bed$", "\\1", basename(beds))
ok_tx <- tx_of %in% txs; beds <- beds[ok_tx]; tx_of <- tx_of[ok_tx]
chr_of <- sub("_.*", "", vapply(beds, function(f) fread(f, nrows = 1, header = FALSE)$V2, ""))
cat(sprintf("transcripts with CpG bed + expression in both: %d\n", length(beds)))
chr_rows <- split(seq_along(gr), as.character(seqnames(gr)))

fit_one <- function(bed, tx, pos2i, Xb, Xa) {
    ids <- fread(bed, header = FALSE)$V2
    li <- pos2i[sub(".*_", "", ids)]; li <- li[!is.na(li)]
    if (length(li) < 200) return(NULL)
    yb <- c(trB[tx, s_tr], teB[tx, s_te]); ya <- c(trA[tx, s_tr], teA[tx, s_te])
    if (!all(is.finite(yb)) || !all(is.finite(ya))) return(NULL)
    xb <- Xb[li, , drop = FALSE]; xa <- Xa[li, , drop = FALSE]
    out <- lapply(LAMBDA, function(l) {
        X <- l * xa + (1 - l) * xb; y <- l * ya + (1 - l) * yb
        Xtr <- t(X[, jtr, drop = FALSE]); Xte <- t(X[, jte, drop = FALSE])
        ok <- apply(Xtr, 2, function(z) all(is.finite(z)) && sd(z) > 0) &
              apply(Xte, 2, function(z) all(is.finite(z)))
        Xtr <- Xtr[, ok, drop = FALSE]; Xte <- Xte[, ok, drop = FALSE]
        if (ncol(Xtr) < K) return(NULL)
        ytr <- y[jtr]; yte <- y[jte]
        sel <- order(abs(cor(Xtr, ytr)), decreasing = TRUE)[1:K]
        f <- lm(ytr ~ ., data = as.data.frame(Xtr[, sel, drop = FALSE]))
        res <- as.numeric(yte - predict(f, newdata = as.data.frame(Xte[, sel, drop = FALSE])))
        list(res = data.table(sample = s_te, tx = tx, lambda = l, res = res),
             r2  = data.table(tx = tx, lambda = l, test_r2 = 1 - sum(res^2) / sum((yte - mean(yte))^2)))
    })
    if (any(vapply(out, is.null, TRUE))) return(NULL)
    list(res = rbindlist(lapply(out, `[[`, "res")), r2 = rbindlist(lapply(out, `[[`, "r2")))
}

ckdir <- file.path(outdir, "per_chr"); dir.create(ckdir, showWarnings = FALSE)
for (ch in names(chr_rows)) {
    sel_tx <- which(chr_of == ch); if (!length(sel_tx)) next
    fres <- file.path(ckdir, paste0(ch, "_res.tsv.gz")); fr2 <- file.path(ckdir, paste0(ch, "_r2.tsv.gz"))
    if (file.exists(fres) && file.exists(fr2)) { cat(sprintf("%-6s already done, skipping\n", ch)); next }
    ri <- chr_rows[[ch]]; t0 <- Sys.time()
    Xb <- as.matrix(MB[ri, iA, drop = FALSE]); Xa <- as.matrix(MA[ri, iB, drop = FALSE])
    pos2i <- setNames(seq_along(ri), as.character(start(gr)[ri]))
    o <- Filter(Negate(is.null), mclapply(sel_tx, function(k) fit_one(beds[k], tx_of[k], pos2i, Xb, Xa),
                                          mc.cores = ncore))
    fwrite(rbindlist(lapply(o, `[[`, "r2")), paste0(fr2, ".tmp"), sep = "\t", compress = "gzip")
    fwrite(rbindlist(lapply(o, `[[`, "res")), paste0(fres, ".tmp"), sep = "\t", compress = "gzip")
    file.rename(paste0(fr2, ".tmp"), fr2); file.rename(paste0(fres, ".tmp"), fres)   # atomic: no half-written checkpoints
    rm(Xb, Xa); gc()
    cat(sprintf("%-6s %5d tx -> %5d modelled  [%s]\n", ch, length(sel_tx), length(o),
                format(round(difftime(Sys.time(), t0, units = "mins"), 1)))); flush.console()
}
res <- rbindlist(lapply(list.files(ckdir, "_res\\.tsv\\.gz$", full.names = TRUE), fread))
r2  <- rbindlist(lapply(list.files(ckdir, "_r2\\.tsv\\.gz$",  full.names = TRUE), fread))
fwrite(res, file.path(outdir, "lambda_resid_long.tsv.gz"), sep = "\t")
fwrite(r2,  file.path(outdir, "lambda_test_r2.tsv.gz"),   sep = "\t")

pn <- read.delim("resources/purity/pames_purity_NT.tsv"); pm <- read.delim("resources/purity/pames_purity_M.tsv")
pur <- c(setNames(pn$Purity[pn$Condition == "Tumour"], pn$Sample[pn$Condition == "Tumour"]), setNames(pm[[2]], pm[[1]]))
pp  <- read.delim("diagnostics/purple_purity_by_met_site.tsv"); purple <- setNames(pp$purity, pp$sample)
agg <- res[, .(absres = mean(abs(res))), by = .(sample, lambda)]
agg[, `:=`(pames = pur[sample], purple = purple[sample],
           group = fifelse(grepl("^DTB", sample), "M", fifelse(grepl("^T", sample), "T", "N")))]
fwrite(agg, file.path(outdir, "lambda_resid_per_sample.tsv"), sep = "\t")

cat(sprintf("\ntranscripts modelled: %d\n\n", uniqueN(r2$tx)))
cat("lambda | median test R2 | |res| vs PAMES: all     T      M  | vs PURPLE: M\n")
for (l in LAMBDA) {
    a <- agg[lambda == l & !is.na(pames) & group != "N"]
    cc <- function(s, x) if (nrow(s) > 4) cor(s$absres, s[[x]], use = "complete.obs") else NA
    cat(sprintf("%5.2f  |     %.4f     |        %+.3f  %+.3f  %+.3f |         %+.3f\n", l,
        median(r2[lambda == l, test_r2]), cc(a, "pames"), cc(a[group == "T"], "pames"),
        cc(a[group == "M"], "pames"), cc(a[group == "M"], "purple")))
}
