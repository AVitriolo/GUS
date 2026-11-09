# =============================================================================
# purify_tumour_NT.R
# CAMDAC deconvolution of CPGEA Tumour methylation (Larose Cadieux et al.,
# bioRxiv 2020.11.03.366252, eq. 2), genome-wide:
#
#   b_t = [b_bulk*(rho*n_t + (1-rho)*n_n) - (1-rho)*n_n*b_n] / (rho*n_t)
#
#   rho : tumour purity
#   n_t : tumour copy number - 2 * 2^log2_ratio inside a GISTIC focal peak
#         called in that sample (cna_peaks_NT.tsv), 2 everywhere else
#   n_n : 2
#   b_n : the patient-matched CPGEA Normal
#
# Every Tumour column is replaced by its purified beta b_t (clipped to [0,1])
# and M is recomputed from it. Normal and metastasis columns are copied as is.
# Reads the integrated h5 block by block and writes a NEW h5 of the same shape.
# =============================================================================
suppressPackageStartupMessages({ library(GenomicRanges); library(rhdf5) })

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

in_h5          <- args$in_h5
out_h5         <- args$out_h5
rowranges_path <- args$rowranges
coldata_path   <- args$coldata
peaks_path     <- args$peaks
purity_path    <- args$purity
min_cn         <- as.numeric(args$min_cn)
block_size     <- as.integer(if (is.null(args$block_size)) 50000 else args$block_size)
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
rho   <- setNames(rho[[grep("^purity$", colnames(rho), ignore.case = TRUE)]],
                  rownames(rho))                 # PUREE: "purity", PAMES: "Purity"

n_rows <- length(rr); n_cols <- nrow(cd)
dims   <- h5ls(in_h5)
stopifnot(identical(dims[dims$name == "beta", "dim"], sprintf("%d x %d", n_rows, n_cols)))

# ---- every Tumour: column, matched Normal column, purity ----
tumours <- rownames(cd)[cd$condition == "Tumour"]
normals <- rownames(cd)[cd$condition == "Normal"]

if (length(setdiff(tumours, names(rho))))
    stop("No purity for: ", paste(setdiff(tumours, names(rho)), collapse = ", "))
t_col <- match(tumours, rownames(cd))
n_col <- match(normals[match(cd[tumours, "patient"], cd[normals, "patient"])], rownames(cd))
if (anyNA(n_col))
    stop("No matched Normal for: ", paste(tumours[is.na(n_col)], collapse = ", "))
t_rho <- rho[tumours]

# ---- per-CpG tumour copy number: rows inside a called peak -> 2*2^log2 ----
peaks     <- peaks[peaks$sample %in% tumours, ]
peaks$n_t <- pmax(N_N * 2^peaks$log2_ratio, min_cn)
peak_gr   <- GRanges(peaks$seqnames, IRanges(peaks$start, peaks$end))
hits      <- findOverlaps(rr, peak_gr)
cn_row    <- queryHits(hits)
cn_tum    <- match(peaks$sample[subjectHits(hits)], tumours)
cn_val    <- peaks$n_t[subjectHits(hits)]
if (anyDuplicated(cbind(cn_row, cn_tum)))
    stop("A CpG falls in two peaks called in the same sample.")

cat(sprintf("%d Tumours (%d with >=1 peak call), %d CpG x sample copy-number entries\n",
            length(tumours), length(unique(peaks$sample)), length(cn_row)))

# ---- output h5, same layout as dataset_integration_split.R ----
if (file.exists(out_h5)) file.remove(out_h5)
h5createFile(out_h5)
for (nm in c("beta", "M"))
    h5createDataset(out_h5, nm, dims = c(n_rows, n_cols),
                    chunk = c(1000, n_cols), level = 0, storage.mode = "double")

n_t_all   <- numeric(length(tumours))            # per-sample running sums for the summary
n_clip    <- numeric(length(tumours))
sum_delta <- numeric(length(tumours))

starts <- seq(1L, n_rows, by = block_size)
for (b in seq_along(starts)) {
    rows <- seq(starts[b], min(starts[b] + block_size - 1L, n_rows))
    beta <- h5read(in_h5, "beta", index = list(rows, NULL))
    M    <- h5read(in_h5, "M",    index = list(rows, NULL))

    n_t <- matrix(N_N, length(rows), length(tumours))
    sel <- cn_row >= rows[1] & cn_row <= rows[length(rows)]
    n_t[cbind(cn_row[sel] - rows[1] + 1L, cn_tum[sel])] <- cn_val[sel]

    b_bulk <- beta[, t_col, drop = FALSE]
    b_n    <- beta[, n_col, drop = FALSE]
    r      <- matrix(t_rho, length(rows), length(tumours), byrow = TRUE)
    b_t    <- (b_bulk * (r * n_t + (1 - r) * N_N) - (1 - r) * N_N * b_n) / (r * n_t)

    n_clip    <- n_clip + colSums(b_t < 0 | b_t > 1)
    b_t       <- pmin(pmax(b_t, 0), 1)
    sum_delta <- sum_delta + colSums(abs(b_t - b_bulk))
    n_t_all   <- n_t_all + length(rows)

    beta[, t_col] <- b_t
    M[, t_col]    <- beta_to_m(b_t)
    h5write(beta, out_h5, "beta", index = list(rows, NULL))
    h5write(M,    out_h5, "M",    index = list(rows, NULL))

    if (b %% 20 == 0 || b == length(starts))
        cat(sprintf("[%s] block %d/%d (rows up to %d)\n",
                    format(Sys.time(), "%H:%M:%S"), b, length(starts), rows[length(rows)]))
}
h5closeAll()

summary <- data.frame(sample = tumours, normal = rownames(cd)[n_col], rho = t_rho,
                      frac_clipped = n_clip / n_t_all,
                      mean_abs_delta = sum_delta / n_t_all)
write.table(summary, output_log, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Done: median fraction clipped %.3f, median |purified - bulk| %.4f\n",
            median(summary$frac_clipped), median(summary$mean_abs_delta)))
