# =============================================================================
# camdac_correct_NT.R
# Copy-number correction of CPGEA Tumour beta values at the GISTIC focal
# peaks, using the CAMDAC mixture model (Larose Cadieux et al., bioRxiv
# 2020.11.03.366252):
#
#   b_bulk = [rho*n_t*b_t + (1-rho)*n_n*b_n] / [rho*n_t + (1-rho)*n_n]
#   b_t    = [b_bulk*(rho*n_t + (1-rho)*n_n) - (1-rho)*n_n*b_n] / (rho*n_t)
#
#   rho : tumour purity (PUREE)
#   n_t : tumour copy number in the peak, 2 * 2^log2_ratio (the CPGEA CN
#         values are purity/ploidy-adjusted: they sit on log2(CN/2) steps)
#   n_n : 2
#   b_n : the patient-matched CPGEA Normal
#
# Only CpGs inside a peak called in that sample are touched; every other value
# (Normal, metastasis, copy-neutral Tumour CpGs) is left as the bulk beta.
# To keep corrected CpGs on the same scale as their uncorrected neighbours,
# the output is the bulk beta the sample WOULD have if the region were
# diploid (n_t = 2), i.e. rho*b_t + (1-rho)*b_n, not the purified b_t.
#
# Works IN PLACE on a copy of the integrated h5 (made by the rule): only the
# rows inside peaks are read and rewritten, for both the beta and M assays.
# =============================================================================
suppressPackageStartupMessages({ library(GenomicRanges); library(rhdf5) })

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

h5_path        <- args$h5                 # copy of integrated_wgbs.h5, modified in place
rowranges_path <- args$rowranges
coldata_path   <- args$coldata
peaks_path     <- args$peaks
purity_path    <- args$purity
min_cn         <- as.numeric(args$min_cn)
output_log     <- args$output_summary

N_N <- 2

beta_to_m <- function(beta, eps = 1e-6) {      # same as dataset_integration_split.R
    beta <- pmin(pmax(beta, eps), 1 - eps)
    log2(beta / (1 - beta))
}

rr    <- readRDS(rowranges_path)
cd    <- as.data.frame(readRDS(coldata_path))
peaks <- read.delim(peaks_path)
rho   <- read.delim(purity_path, row.names = 1)
rho   <- data.frame(purity = rho[[grep("^purity$", colnames(rho), ignore.case = TRUE)]],
                    row.names = rownames(rho))     # PUREE: "purity", PAMES: "Purity"

dims <- h5ls(h5_path)
dims <- dims[dims$name == "beta", "dim"]
stopifnot(identical(dims, sprintf("%d x %d", length(rr), nrow(cd))))
if (is.unsorted(order(seqnames(rr), start(rr))))
    stop("rowRanges are not sorted: peaks would not map to contiguous rows.")

# ---- Tumour samples with calls, their purity and matched Normal ----
tumours <- rownames(cd)[cd$condition == "Tumour"]
normals <- rownames(cd)[cd$condition == "Normal"]
peaks   <- peaks[peaks$sample %in% tumours, ]

missing_rho <- setdiff(unique(peaks$sample), rownames(rho))
if (length(missing_rho)) stop("No purity for: ", paste(missing_rho, collapse = ", "))

normal_of <- setNames(normals[match(cd[unique(peaks$sample), "patient"],
                                    cd[normals, "patient"])],
                      unique(peaks$sample))
if (anyNA(normal_of))
    stop("No matched Normal for: ", paste(names(normal_of)[is.na(normal_of)], collapse = ", "))

peaks$col   <- match(peaks$sample, rownames(cd))
peaks$ncol  <- match(normal_of[peaks$sample], rownames(cd))
peaks$rho   <- rho[peaks$sample, "purity"]
peaks$n_t   <- pmax(N_N * 2^peaks$log2_ratio, min_cn)

cat(sprintf("%d calls in %d Tumour samples (of %d); n_t range %.2f-%.2f\n",
            nrow(peaks), length(unique(peaks$sample)), length(tumours),
            min(peaks$n_t), max(peaks$n_t)))

# ---- correct peak by peak (each peak = one contiguous block of rows) ----
peak_gr <- unique(GRanges(peaks$seqnames, IRanges(peaks$start, peaks$end), peak = peaks$peak))
summary <- list()

for (i in seq_along(peak_gr)) {
    pk   <- peak_gr[i]
    rows <- which(overlapsAny(rr, pk))
    if (!length(rows)) next
    stopifnot(all(diff(rows) == 1L))

    calls <- peaks[peaks$peak == pk$peak, ]
    idx   <- list(rows, sort(unique(c(calls$col, calls$ncol))))
    beta  <- h5read(h5_path, "beta", index = idx)
    cols  <- idx[[2]]

    for (k in seq_len(nrow(calls))) {
        c_ <- calls[k, ]
        b_bulk <- beta[, cols == c_$col]
        b_n    <- beta[, cols == c_$ncol]
        mix    <- c_$rho * c_$n_t + (1 - c_$rho) * N_N
        b_t    <- (b_bulk * mix - (1 - c_$rho) * N_N * b_n) / (c_$rho * c_$n_t)
        b_t    <- pmin(pmax(b_t, 0), 1)
        b_new  <- c_$rho * b_t + (1 - c_$rho) * b_n

        summary[[length(summary) + 1]] <- data.frame(
            sample = c_$sample, peak = c_$peak, n_cpgs = length(rows),
            rho = c_$rho, n_t = c_$n_t,
            frac_clipped   = mean(b_t %in% c(0, 1)),
            mean_abs_delta = mean(abs(b_new - b_bulk)))

        # write the Tumour column only; its Normal is never modified
        h5write(matrix(b_new, ncol = 1),            h5_path, "beta", index = list(rows, c_$col))
        h5write(matrix(beta_to_m(b_new), ncol = 1), h5_path, "M",    index = list(rows, c_$col))
    }
    cat(sprintf("[%s] %s: %d CpGs x %d samples corrected\n",
                pk$peak, as.character(pk), length(rows), nrow(calls)))
}
h5closeAll()

summary <- do.call(rbind, summary)
write.table(summary, output_log, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Done: %d (sample, peak) corrections; median |delta beta| = %.4f\n",
            nrow(summary), median(summary$mean_abs_delta)))
