# =============================================================================
# run_bmind.R
# Reference-free tumour RNA purification with bMIND (MIND package, Wang et al.
# Nat Commun 2021): bulk log2(CPM + 1) + per-sample fractions
# (Tumour = purity, Other = 1 - purity), profile = covariance = NULL, so the
# prior is learned from the samples themselves. Used for CPGEA and for the
# metastases (one fit per biopsy site), with purity from methylation.
#
#   --fit_samples : samples the model is fitted on (comma list or file)
#   --keep        : samples written out (default: all fitted). For a test set,
#                   fit train + test together and keep only test, so test
#                   borrows the train fit and nothing flows back to train.
#   CPM is computed over ALL transcripts of each sample before subsetting.
#
# Output: <out_prefix>_tumour_log2cpm.tsv.gz (bMIND tumour estimate, A[, "Tumour", ]),
#         <out_prefix>_tumour_se.tsv.gz, <out_prefix>_bulk_log2cpm.tsv.gz
# =============================================================================
suppressPackageStartupMessages({ library(MIND); library(data.table) })

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
counts_path <- args$counts
tx_path     <- args$transcripts
purity_path <- args$purity
out_prefix  <- args$out_prefix
ncore       <- as.integer(if (is.null(args$ncore)) 16 else args$ncore)
read_ids    <- function(x) if (file.exists(x)) readLines(x) else strsplit(x, ",")[[1]]
fit_samples <- read_ids(args$fit_samples)
keep        <- if (is.null(args$keep)) fit_samples else read_ids(args$keep)
stopifnot(all(keep %in% fit_samples))
dir.create(dirname(out_prefix), showWarnings = FALSE, recursive = TRUE)

pm <- read.delim(purity_path, row.names = 1)
pi <- setNames(pm[[grep("^purity$", colnames(pm), ignore.case = TRUE)]], rownames(pm))
if (length(setdiff(fit_samples, names(pi))))
    stop("No purity for: ", paste(setdiff(fit_samples, names(pi)), collapse = ", "))
pi <- pi[fit_samples]

raw <- fread(counts_path)
if (length(setdiff(fit_samples, names(raw)))) stop("Samples missing from counts: ",
    paste(setdiff(fit_samples, names(raw)), collapse = ", "))
m <- as.matrix(raw[, ..fit_samples]); rownames(m) <- raw[[1]]
lcpm <- log2(t(t(m) / colSums(m)) * 1e6 + 1)

tx <- intersect(fread(tx_path, select = 1)[[1]], rownames(lcpm))
X  <- lcpm[tx, , drop = FALSE]
X  <- X[apply(X, 1, sd) > 0, , drop = FALSE]          # bMIND needs variation across samples
stopifnot(max(X) < 50)                                # already log2: bMIND must not re-transform
cat(sprintf("fit on %d samples (keep %d) | %d transcripts | purity median %.2f (%.2f-%.2f)\n",
            ncol(X), length(keep), nrow(X), median(pi), min(pi), max(pi)))

frac <- data.frame(Tumour = pi, Other = 1 - pi, row.names = fit_samples)
t0  <- Sys.time()
fit <- bMIND(X, frac = frac, ncore = ncore)
cat(sprintf("bMIND: %d of %d transcripts returned in %.1f min\n", dim(fit$A)[1], nrow(X),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

A  <- fit$A[, "Tumour", keep, drop = FALSE][, 1, ]
SE <- fit$SE[, "Tumour", keep, drop = FALSE][, 1, ]
if (length(keep) == 1) { A <- matrix(A, ncol = 1); SE <- matrix(SE, ncol = 1) }
dimnames(A) <- dimnames(SE) <- list(dimnames(fit$A)[[1]], keep)

fwrite(data.table(transcript_id = rownames(A), A), paste0(out_prefix, "_tumour_log2cpm.tsv.gz"), sep = "\t")
fwrite(data.table(transcript_id = rownames(SE), SE), paste0(out_prefix, "_tumour_se.tsv.gz"), sep = "\t")
fwrite(data.table(transcript_id = rownames(A), X[rownames(A), keep, drop = FALSE]),
       paste0(out_prefix, "_bulk_log2cpm.tsv.gz"), sep = "\t")
fwrite(data.table(sample = keep, purity = pi[keep]), paste0(out_prefix, "_purity.tsv"), sep = "\t")

shift <- A - X[rownames(A), keep, drop = FALSE]
cat(sprintf("written %d samples | median |shift| %.3f | transcripts with |median shift| > 1: %d\n",
            length(keep), median(abs(shift)), sum(abs(apply(shift, 1, median)) > 1)))
