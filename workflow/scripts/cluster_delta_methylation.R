# For each cluster in the selected-features BED, compute mean methylation
# in two conditions and their delta, from the integrated (combat) H5.
#
# Conditions (--cond_1 / --cond_2) accept either the short code or the
# colData$condition value:
#   N / Normal      T / Tumour      D / metastasis
# For a metastasis condition, --site_1 / --site_2 optionally restrict to a
# metastatic site (colData$metastatis_site, e.g. Bone / Liver / Lymph_node).
#
# Sample sizes are matched by random subsampling to the smaller group, so the
# delta is not driven by unequal n. set.seed(123) keeps it reproducible.
#
# Input BED (one row per CpG; from selected_features_to_bed.R --mode=percpg,
# 5-col variant):   chr  start  end  CpGID(chr_pos)  clusterID
#
# Output TSV (one row per cluster):
#   clusterID  chr  start  end  n_CpGs  <label1>  <label2>  delta
#   delta = mean(cond_2) - mean(cond_1), in M-value space (H5 is M-only).

set.seed(123)
suppressPackageStartupMessages({ library(dplyr) })
options(scipen = 999)

source("workflow/scripts/helpers/load_h5_rse_integrated.R")

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
input_h5  <- args$input_h5
input_bed <- args$input_bed
cond_1    <- args$cond_1
cond_2    <- args$cond_2
site_1    <- args$site_1          # optional, only for a metastasis condition
site_2    <- args$site_2          # optional
output    <- args$output

stopifnot(!is.null(cond_1), !is.null(cond_2), !is.null(input_bed), !is.null(output))
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)

# map short codes to colData$condition values
cond_map <- c(N = "Normal", T = "Tumour", D = "metastasis",
              Normal = "Normal", Tumour = "Tumour", metastasis = "metastasis")
resolve_cond <- function(x) {
    v <- cond_map[[x]]
    if (is.null(v)) stop(sprintf("unknown condition '%s' (use N/T/D)", x))
    v
}
cond_1_val <- resolve_cond(cond_1)
cond_2_val <- resolve_cond(cond_2)

# ---- load integrated methylation (M assay; CpGs x samples in R's view) ----
obj  <- load_h5_rse_integrated(input_h5)
M    <- obj$M
rownames(M) <- obj$CpGs$CpGID
cd   <- SummarizedExperiment::colData(obj$rse)
colnames(M) <- rownames(cd)

cat(sprintf("M: %d CpGs x %d samples\n", nrow(M), ncol(M)))

# ---- pick samples for a condition (+ optional site) ----
pick_samples <- function(cond_val, site) {
    keep <- cd$condition == cond_val
    if (!is.null(site)) {
        if (cond_val != "metastasis")
            stop(sprintf("--site given but condition '%s' is not metastasis", cond_val))
        keep <- keep & !is.na(cd$metastatis_site) & cd$metastatis_site == site
    }
    rownames(cd)[keep]
}

s1 <- pick_samples(cond_1_val, site_1)
s2 <- pick_samples(cond_2_val, site_2)

label_1 <- if (!is.null(site_1)) paste0(cond_1, "_", site_1) else cond_1
label_2 <- if (!is.null(site_2)) paste0(cond_2, "_", site_2) else cond_2

cat(sprintf("cond_1 = %s : %d samples\n", label_1, length(s1)))
cat(sprintf("cond_2 = %s : %d samples\n", label_2, length(s2)))
if (length(s1) == 0 || length(s2) == 0) stop("one condition has no samples")

# ---- size-match by subsampling the larger group ----
n <- min(length(s1), length(s2))
if (length(s1) > n) s1 <- sample(s1, n)
if (length(s2) > n) s2 <- sample(s2, n)
cat(sprintf("size-matched to %d samples per group\n", n))

# ---- load the selected-features BED (one row per CpG) ----
bed <- read.table(input_bed, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(bed) <- c("chr", "start", "end", "CpGID", "clusterID")
cat(sprintf("BED: %d CpG rows, %d clusters\n",
            nrow(bed), length(unique(bed$clusterID))))
# ---- unique CpGs actually needed ----
bed <- bed[bed$CpGID %in% rownames(M), , drop = FALSE]
if (nrow(bed) == 0) stop("no BED CpGs found in the H5 (CpGID mismatch)")

uniq_cpgs <- unique(bed$CpGID)
cat(sprintf("unique CpGs to read: %d\n", length(uniq_cpgs)))

# realize just those rows for the two sample groups (as ordinary matrices)
row_idx <- match(uniq_cpgs, rownames(M))
M1 <- as.matrix(M[row_idx, s1, drop = FALSE])
M2 <- as.matrix(M[row_idx, s2, drop = FALSE])

mean_1 <- rowMeans(M1, na.rm = TRUE)
mean_2 <- rowMeans(M2, na.rm = TRUE)
names(mean_1) <- uniq_cpgs
names(mean_2) <- uniq_cpgs

# map per-CpG means back onto every BED row
bed$m1 <- mean_1[bed$CpGID]
bed$m2 <- mean_2[bed$CpGID]

# ---- aggregate per cluster ----
clusters <- bed %>%
    group_by(clusterID) %>%
    summarise(
        chr    = first(chr),
        start  = min(start),
        end    = max(end),
        n_CpGs = n(),
        cond_1 = mean(m1, na.rm = TRUE),
        cond_2 = mean(m2, na.rm = TRUE),
        delta  = mean(m2, na.rm = TRUE) - mean(m1, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    rename(!!label_1 := cond_1, !!label_2 := cond_2) %>%
    arrange(desc(abs(delta)))

write.table(clusters, file = output, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("written %d clusters -> %s\n", nrow(clusters), output))
cat(sprintf("delta = %s minus %s (M-value space)\n", label_2, label_1))