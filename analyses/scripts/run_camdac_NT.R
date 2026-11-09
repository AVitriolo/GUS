# =============================================================================
# run_camdac_NT.R
# Purified CPGEA Tumour methylation with the CAMDAC package (v0.2.0,
# VanLoo-lab/CAMDAC), starting from per-CpG read counts instead of BAMs:
#
#   per Tumour / matched-Normal pair
#     meth   : <id>.m.csv.gz built from the CPGEA source SE (beta, cov):
#              chrom, start, end, m, cov, M, UM, BAF, cg_snp
#              (no SNP information -> BAF / cg_snp NA, so CAMDAC uses total CN)
#     cna    : genome-wide segments chrom, start, end, major_cn, minor_cn,
#              purity, ploidy. Called GISTIC peaks carry the sample's copy
#              number (2 * 2^log2_ratio, total CN on major_cn); every other
#              base is CN 2 (1 + 1). CpGs outside all segments would be dropped
#              by CAMDAC, hence the tiling.
#     attach_output(meth, cna) -> CAMDAC::cmain_deconvolve_methylation()
#
#   chrX / chrY are set to CN 2 like the autosomes: CAMDAC assumes a normal CN
#   of 2 everywhere, so in these male samples the tumour/normal CN ratio on X/Y
#   (1:1) is only preserved by giving the tumour 2 as well.
# =============================================================================
suppressPackageStartupMessages({
    library(CAMDAC); library(SummarizedExperiment); library(data.table); library(GenomicRanges)
})

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

source_se    <- args$source_se           # resources/cpgea_wgbs_with_coverage_hg38/se.rds
coldata_path <- args$coldata             # integrated colData (defines the T/N pairs used)
peaks_path   <- args$peaks
purity_path  <- args$purity
chrom_sizes  <- args$chrom_sizes
outdir       <- args$outdir
n_cores      <- as.integer(if (is.null(args$n_cores)) 1 else args$n_cores)
only         <- if (is.null(args$tumours)) NULL else strsplit(args$tumours, ",")[[1]]

CHROMS <- paste0("chr", c(1:22, "X", "Y"))

# ---- pairs, purity, peaks ----
cd      <- as.data.frame(readRDS(coldata_path))
tumours <- rownames(cd)[cd$condition == "Tumour"]
normals <- rownames(cd)[cd$condition == "Normal"]
pairs   <- data.frame(tumour = tumours,
                      normal = normals[match(cd[tumours, "patient"], cd[normals, "patient"])])
if (anyNA(pairs$normal)) stop("Tumours without a matched Normal.")
if (!is.null(only)) pairs <- pairs[pairs$tumour %in% only, ]

rho <- read.delim(purity_path, row.names = 1)
rho <- setNames(rho[[grep("^purity$", colnames(rho), ignore.case = TRUE)]], rownames(rho))
if (length(setdiff(pairs$tumour, names(rho)))) stop("Tumours without purity.")

peaks <- read.delim(peaks_path)
sizes <- fread(chrom_sizes, header = FALSE, col.names = c("chrom", "size"))
sizes <- sizes[chrom %in% CHROMS]

# ---- source counts: loaded as an HDF5 SummarizedExperiment, like
#      dataset_integration_split.R. Do NOT read assays.h5 directly: the SE's
#      assays are a row-reordered DelayedSubset of it (rows >= 2,375,160). ----
se <- HDF5Array::loadHDF5SummarizedExperiment(dirname(source_se))
gr <- rowRanges(se)
pos <- data.table(chrom = as.character(seqnames(gr)), start = start(gr))
keep_chr <- pos$chrom %in% CHROMS

config <- CamConfig(outdir = outdir, bsseq = "wgbs", lib = "pe", build = "hg38",
                    n_cores = n_cores, cna_caller = "ascat")
stage  <- file.path(outdir, "_stage")
dir.create(stage, showWarnings = FALSE, recursive = TRUE)

write_meth <- function(id, path) {
    j    <- match(id, colnames(se))
    beta <- as.vector(assay(se, "beta")[, j])
    cov  <- as.vector(assay(se, "cov")[, j])
    ok   <- keep_chr & !is.na(beta) & !is.na(cov) & cov > 0
    M    <- round(beta[ok] * cov[ok])
    dt   <- data.table(chrom = pos$chrom[ok], start = pos$start[ok], end = pos$start[ok] + 1L,
                       m = M / cov[ok], cov = cov[ok], M = M, UM = cov[ok] - M,
                       BAF = NA_real_, cg_snp = NA_character_)
    fwrite(dt, path)
    nrow(dt)
}

write_cna <- function(tumour, path) {
    p <- peaks[peaks$sample == tumour, ]
    p <- p[order(match(p$seqnames, CHROMS), p$start), ]
    segs <- rbindlist(lapply(seq_len(nrow(sizes)), function(i) {
        ch <- sizes$chrom[i]; pk <- p[p$seqnames == ch, ]
        cuts <- c(1, rbind(pk$start, pk$end + 1), sizes$size[i] + 1)
        s <- data.table(chrom = ch, start = head(cuts, -1), end = tail(cuts, -1) - 1,
                        major_cn = 1, minor_cn = 1)
        if (nrow(pk)) {                       # every 2nd segment is a called peak
            idx <- seq(2, by = 2, length.out = nrow(pk))
            s$major_cn[idx] <- 2 * 2^pk$log2_ratio
            s$minor_cn[idx] <- 0
        }
        s[end >= start]
    }))
    segs[, `:=`(purity = rho[[tumour]], ploidy = 2)]
    fwrite(segs, path, sep = "\t")
    nrow(segs)
}

for (k in seq_len(nrow(pairs))) {
    t_id <- pairs$tumour[k]; n_id <- pairs$normal[k]
    t0   <- Sys.time()
    tumour <- CamSample(id = t_id, sex = "XY", patient_id = t_id)
    normal <- CamSample(id = n_id, sex = "XY", patient_id = t_id)

    if (file.exists(CAMDAC:::get_fpath(tumour, config, "pure"))) {
        cat(sprintf("[%s] pure exists, skipping\n", t_id)); next
    }

    f_t <- file.path(stage, paste0(t_id, ".m.csv.gz"))
    f_n <- file.path(stage, paste0(n_id, ".m.csv.gz"))
    f_c <- file.path(stage, paste0(t_id, ".cna.txt"))
    n_t <- write_meth(t_id, f_t); n_n <- write_meth(n_id, f_n); n_s <- write_cna(t_id, f_c)

    attach_output(tumour, config, "meth", f_t)
    attach_output(normal, config, "meth", f_n)
    attach_output(tumour, config, "cna",  f_c)
    file.remove(f_t, f_n, f_c)

    cmain_deconvolve_methylation(tumour, normal, config)

    # the attached bulk meth copies are CAMDAC inputs only; keep just the output
    file.remove(CAMDAC:::get_fpath(tumour, config, "meth"), CAMDAC:::get_fpath(normal, config, "meth"))
    cat(sprintf("[%s/%s] %d/%d CpGs, %d CNA segments, rho %.2f -> %s (%.1f min)\n",
                t_id, n_id, n_t, n_n, n_s, rho[[t_id]],
                CAMDAC:::get_fpath(tumour, config, "pure"),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
