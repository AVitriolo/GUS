# Subset of the purified integrated h5 (/home/ieo7293/Prostate_GAS_purified, CAMDAC-purified tumours/metastases,
# bulk normals, v1 CpG universe) to the CpGs of a MethylDriver per-CpG BED, with the purified BETA stored as dataset
# "M" so cluster_delta_methylation.R (which reads "M") works on the beta scale, like the bulk run on
# integrated_wgbs_combat_beta.h5. Writes <out>.h5, <out>_rowRanges.rds, <out>_colData.rds
# Usage: Rscript subset_purified_h5_beta.R <percpg.bed> <out.h5>
suppressPackageStartupMessages({library(HDF5Array); library(data.table); library(GenomicRanges)})
a <- commandArgs(TRUE); bed <- fread(a[1], header = FALSE); out <- a[2]
P <- "/home/ieo7293/Prostate_GAS_purified/resources/integrated_wgbs_hdf5/integrated_wgbs"
rr <- readRDS(paste0(P, "_rowRanges.rds")); cd <- readRDS(paste0(P, "_colData.rds"))
id <- paste0(seqnames(rr), "_", start(rr))
want <- unique(bed$V4); i <- match(want, id); keep <- sort(unique(i[!is.na(i)]))
cat(sprintf("BED CpGs %d | in the purified (v1) universe %d (%.1f%%)\n", length(want), length(keep), 100 * length(keep) / length(want)))
B <- HDF5Array(paste0(P, ".h5"), "beta")          # CpGs x samples
X <- as.matrix(B[keep, , drop = FALSE])
cat(sprintf("extracted %d x %d, NA %d\n", nrow(X), ncol(X), sum(is.na(X))))
if (file.exists(out)) file.remove(out)
writeHDF5Array(X, filepath = out, name = "M", chunkdim = c(min(1000, nrow(X)), ncol(X)))
saveRDS(rr[keep], sub("\\.h5$", "_rowRanges.rds", out)); saveRDS(cd, sub("\\.h5$", "_colData.rds", out))
print(table(as.data.frame(cd)$condition))
