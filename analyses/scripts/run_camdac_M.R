# =============================================================================
# run_camdac_M.R
# Purified WCDT metastatic methylation with the CAMDAC package (v0.2.0), the
# metastatic counterpart of run_camdac_NT.R:
#
#   tumour meth : MCRPC WGBS counts (resources/mcrpc_wgbs_hg38), liver excluded
#   infiltrates : the patient's own benign-adjacent WGBS if available (9
#                 patients), otherwise the pooled benign normals of the same
#                 biopsy site (counts summed across that site's normals; the
#                 CAMDAC fallback of a tissue-matched normal panel)
#   purity      : PAMES (pames_purity_M.tsv)
#   CNA         : copycat segments converted to tumour copy number. copycat is
#                 relative to the sample average (REF ~2 whatever the ploidy), so
#                 with R = value / 2, purity p and WCDT ploidy psi:
#                   n_t = [R (p psi + 2 (1 - p)) - 2 (1 - p)] / p
#                 for called GAIN / LOSS segments; REF segments and gaps -> psi
#                 (R = 1). Total CN on major_cn (no allele-specific data).
# =============================================================================
suppressPackageStartupMessages({
    library(CAMDAC); library(SummarizedExperiment); library(data.table); library(GenomicRanges)
})

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
mcrpc_path   <- args$mcrpc
coldata_path <- args$coldata
exclude_path <- args$exclude
normal_dir   <- args$normal_dir           # merged benign files (chrom, pos, M, U)
normal_sites <- args$normal_sites         # tsv: normal, site
purity_path  <- args$purity
copycat_dir  <- args$copycat_dir
copycat_tab  <- args$copycat_samples      # tsv: file, biopsy, sample, ploidy
chrom_sizes  <- args$chrom_sizes
outdir       <- args$outdir
min_cn       <- as.numeric(if (is.null(args$min_cn)) 0.1 else args$min_cn)
only         <- if (is.null(args$tumours)) NULL else strsplit(args$tumours, ",")[[1]]

CHROMS <- paste0("chr", c(1:22, "X", "Y"))

cd   <- as.data.frame(readRDS(coldata_path))
mets <- setdiff(rownames(cd)[cd$condition == "metastasis"], readLines(exclude_path))
if (!is.null(only)) mets <- intersect(mets, only)

pm  <- read.delim(purity_path, row.names = 1)
rho <- setNames(pm[[grep("^purity$", colnames(pm), ignore.case = TRUE)]], rownames(pm))
cc  <- fread(copycat_tab)
ns  <- fread(normal_sites)
ns[, patient := gsub("-", "_", sub("^NT-(DTB-[0-9]+).*$", "\\1", normal))]
sizes <- fread(chrom_sizes, header = FALSE, col.names = c("chrom", "size"))[chrom %in% CHROMS]

se <- HDF5Array::loadHDF5SummarizedExperiment(mcrpc_path)
gr <- rowRanges(se)
pos <- data.table(chrom = as.character(seqnames(gr)), start = start(gr))
keep_chr <- pos$chrom %in% CHROMS

config <- CamConfig(outdir = outdir, bsseq = "wgbs", lib = "pe", build = "hg38", cna_caller = "ascat")
stage  <- file.path(outdir, "_stage")
dir.create(stage, showWarnings = FALSE, recursive = TRUE)

as_camdac <- function(chrom, start, M, cov) data.table(
    chrom = chrom, start = start, end = start + 1L, m = M / cov, cov = cov, M = M, UM = cov - M,
    BAF = NA_real_, cg_snp = NA_character_)

write_tumour_meth <- function(id, path) {
    j    <- match(id, colnames(se))
    beta <- as.vector(assay(se, "beta")[, j])
    cov  <- as.vector(assay(se, "cov")[, j])
    ok   <- keep_chr & !is.na(beta) & !is.na(cov) & cov > 0
    M    <- round(beta[ok] * cov[ok])
    fwrite(as_camdac(pos$chrom[ok], pos$start[ok], M, cov[ok]), path)
    sum(ok)
}

# own benign sample, or the pooled normals of the site (built once, cached)
normal_meth <- function(met) {
    own <- ns[patient == met, normal]
    if (length(own)) {
        path <- file.path(stage, paste0(own, ".normal.m.csv.gz"))
        if (!file.exists(path)) {
            x <- fread(file.path(normal_dir, paste0(own, ".merged.tsv.gz")))[M + U > 0]
            fwrite(as_camdac(x$chrom, x$pos, x$M, x$M + x$U), path)
        }
        return(list(id = own, path = path, kind = "own"))
    }
    met_site <- cd[met, "metastatis_site"]          # not `site`: ns has a column of that name
    pool <- ns[site_group == met_site, normal]
    if (!length(pool)) stop("No benign normal for site ", met_site, " (", met, ")")
    id <- paste0("POOL_", met_site); path <- file.path(stage, paste0(id, ".normal.m.csv.gz"))
    if (!file.exists(path)) {
        x <- rbindlist(lapply(pool, function(n) fread(file.path(normal_dir, paste0(n, ".merged.tsv.gz")))))
        x <- x[, .(M = sum(M), U = sum(U)), by = .(chrom, pos)][M + U > 0]
        setorder(x, chrom, pos)
        fwrite(as_camdac(x$chrom, x$pos, x$M, x$M + x$U), path)
    }
    list(id = id, path = path, kind = sprintf("pool of %d", length(pool)))
}

write_cna <- function(met, path) {
    info <- cc[sample == met]
    stopifnot(nrow(info) == 1)
    p <- rho[[met]]; psi <- info$ploidy
    x <- fread(file.path(copycat_dir, info$file), header = FALSE,
               col.names = c("chrom", "start", "end", "status", "value", "strand"))[chrom %in% CHROMS]
    x[, R := ifelse(status == "REF", 1, value / 2)]
    x[, n_t := pmax((R * (p * psi + 2 * (1 - p)) - 2 * (1 - p)) / p, min_cn)]
    segs <- rbindlist(lapply(seq_len(nrow(sizes)), function(i) {
        ch <- sizes$chrom[i]; s <- x[chrom == ch][order(start)]
        # fill gaps (and chromosome ends) with the baseline ploidy
        edges <- c(1, s$end + 1); gaps <- data.table(start = edges, end = c(s$start - 1, sizes$size[i]))
        gaps <- gaps[end >= start][, n_t := psi]
        rbind(s[, .(start, end, n_t)], gaps)[order(start)][, chrom := ch]
    }))
    # copycat segments touch (end b, next start b + 1): a CpG [b, b + 1] would hit
    # both and be duplicated by CAMDAC's overlap join, so each segment after the
    # first starts 1 bp later (no CpG is lost: [b + 1, b + 2] still hits it)
    segs[, start := ifelse(start > 1L, start + 1L, start)]
    segs <- segs[end >= start]
    segs <- segs[, .(chrom, start, end, major_cn = n_t, minor_cn = 0, purity = p, ploidy = psi)]
    fwrite(segs, path, sep = "\t")
    c(n = nrow(segs), called = sum(x$status != "REF"), psi = psi)
}

for (met in mets) {
    t0 <- Sys.time()
    if (!met %in% names(rho)) stop("No purity for ", met)
    tumour <- CamSample(id = met, sex = "XY", patient_id = met)
    if (file.exists(CAMDAC:::get_fpath(tumour, config, "pure"))) { cat(sprintf("[%s] pure exists, skipping\n", met)); next }

    nm <- normal_meth(met)
    normal <- CamSample(id = nm$id, sex = "XY", patient_id = met)
    f_t <- file.path(stage, paste0(met, ".m.csv.gz")); f_c <- file.path(stage, paste0(met, ".cna.txt"))
    n_t <- write_tumour_meth(met, f_t); cn <- write_cna(met, f_c)

    attach_output(tumour, config, "meth", f_t)
    attach_output(normal, config, "meth", nm$path)
    attach_output(tumour, config, "cna",  f_c)
    file.remove(f_t, f_c)

    cmain_deconvolve_methylation(tumour, normal, config)
    file.remove(CAMDAC:::get_fpath(tumour, config, "meth"), CAMDAC:::get_fpath(normal, config, "meth"))
    cat(sprintf("[%s/%s (%s)] %d CpGs, %d segments (%d called), rho %.2f, ploidy %.2f -> %s (%.1f min)\n",
                met, nm$id, nm$kind, n_t, cn[["n"]], cn[["called"]], rho[[met]], cn[["psi"]],
                CAMDAC:::get_fpath(tumour, config, "pure"),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
