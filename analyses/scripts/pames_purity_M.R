# =============================================================================
# pames_purity_M.R
# PAMES purity for the WCDT metastases, mirroring pames_purity.R (same filters,
# parameters and bug-corrected site ranking) with:
#   tumours : MCRPC metastases (resources/mcrpc_wgbs_hg38), minus --exclude
#   controls: the 10 WCDT benign-adjacent samples pooled across sites
#             (bone / lymph node / other soft tissue), strand-merged by
#             merge_benign_M_strands.R
# Informative sites are chosen on all these samples (no train/test split),
# as for CPGEA. Inputs are read chromosome by chromosome.
# =============================================================================
suppressPackageStartupMessages({
    library(HDF5Array); library(SummarizedExperiment); library(GenomicRanges)
    library(data.table); library(PAMES)
})
options(scipen = 999)

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
mcrpc_path      <- args$mcrpc
normal_dir      <- args$normal_dir
coldata_path    <- args$coldata
exclude_path    <- args$exclude
min_cov         <- as.numeric(args$min_cov)
min_sample_frac <- as.numeric(args$min_sample_frac)
min_mean_diff   <- as.numeric(args$min_mean_diff)
max_na_data     <- as.numeric(args$max_na_data)
percentile_low  <- as.numeric(args$percentile_low)
percentile_high <- as.numeric(args$percentile_high)
max_sites       <- as.integer(args$max_sites)
min_distance    <- as.numeric(args$min_distance)
cores           <- as.integer(args$cores)
out_purity      <- args$output_purity
out_sites       <- args$output_info_sites

CHROMS <- paste0("chr", c(1:22, "X", "Y"))

# ---- samples ----
se <- loadHDF5SummarizedExperiment(mcrpc_path)
cd <- as.data.frame(readRDS(coldata_path))
mets <- rownames(cd)[cd$condition == "metastasis"]
if (!is.null(exclude_path)) mets <- setdiff(mets, readLines(exclude_path))
stopifnot(all(mets %in% colnames(se)))
nfiles <- list.files(normal_dir, "\\.merged\\.tsv\\.gz$", full.names = TRUE)
normals <- sub("\\.merged\\.tsv\\.gz$", "", basename(nfiles))
message(sprintf("[%s] Samples: %d metastases, %d benign normals", Sys.time(), length(mets), length(normals)))

nl <- lapply(nfiles, function(f) { x <- fread(f); x[, chrom := factor(chrom, levels = CHROMS)]; x })
names(nl) <- normals

rr   <- rowRanges(se)
chr  <- as.character(seqnames(rr))
cols <- match(mets, colnames(se))

tumor_parts <- list(); control_parts <- list(); ref_parts <- list(); n_cov <- 0L
for (ch in intersect(CHROMS, unique(chr))) {
    idx <- which(chr == ch)
    pos <- start(rr)[idx]
    tcov <- as.matrix(assay(se, "cov")[idx, cols]); tcov[is.na(tcov)] <- 0L
    ncov <- nbeta <- matrix(NA_real_, length(idx), length(normals))
    for (k in seq_along(normals)) {
        x <- nl[[k]][chrom == ch]
        j <- match(pos, x$pos)
        cv <- x$M[j] + x$U[j]
        ncov[, k]  <- ifelse(is.na(cv), 0, cv)
        nbeta[, k] <- x$M[j] / cv
    }
    keep <- which(rowSums(tcov >= min_cov) >= min_sample_frac * length(mets) &
                  rowSums(ncov >= min_cov) >= min_sample_frac * length(normals))
    n_cov <- n_cov + length(keep)
    if (!length(keep)) next
    tbeta <- as.matrix(assay(se, "beta")[idx[keep], cols])
    nb <- nbeta[keep, , drop = FALSE]
    dk <- which(abs(rowMeans(tbeta, na.rm = TRUE) - rowMeans(nb, na.rm = TRUE)) >= min_mean_diff)
    if (length(dk)) {
        key <- paste0(ch, "_", pos[keep][dk])
        tumor_parts[[ch]]   <- `rownames<-`(`colnames<-`(tbeta[dk, , drop = FALSE], mets), key)
        control_parts[[ch]] <- `rownames<-`(`colnames<-`(nb[dk, , drop = FALSE], normals), key)
        ref_parts[[ch]]     <- data.frame(Chromosome = ch, Genomic_Coordinate = pos[keep][dk])
    }
    message(sprintf("[%s] %s: %d CpGs pass coverage, %d pass mean difference", Sys.time(), ch, length(keep), length(dk)))
}
tumor_table   <- do.call(rbind, tumor_parts)
control_table <- do.call(rbind, control_parts)
ref_table     <- do.call(rbind, ref_parts)
message(sprintf("[%s] Candidate CpGs: %d after coverage, %d after mean difference", Sys.time(), n_cov, nrow(tumor_table)))
assertthat::assert_that(nrow(tumor_table) > 0, msg = "No CpG passed the filters.")

# ---- PAMES, as in pames_purity.R ----
auc <- get_AUC(tumor_table, control_table, max_NA_data = max_na_data, cores = cores)
all_sites <- find_informative_sites(tumor_table, control_table, auc, ref_table,
                                    cores = cores, max_sites = max_sites,
                                    min_distance = min_distance,
                                    percentiles = c(percentile_low, percentile_high),
                                    full_info = TRUE)
select_ranked <- function(df, score, n) {
    if (nrow(df) == 0) return(df[0, ])
    df <- df[order(score, decreasing = TRUE), ]
    df$Keep_site <- PAMES:::remove_close_sites(df, n, min_distance)
    head(df[df$Keep_site, ], n)
}
hyper_df <- all_sites[all_sites$Site_type == "hyper", ]
hypo_df  <- all_sites[all_sites$Site_type == "hypo",  ]
hyper_kept <- select_ranked(hyper_df, hyper_df$AUC,    max_sites / 2)
hypo_kept  <- select_ranked(hypo_df,  1 - hypo_df$AUC, max_sites / 2)
info_sites <- rbind(hyper_kept, hypo_kept)[, c("Probe", "Site_type")]
message(sprintf("[%s] Selected %d hyper- and %d hypo-methylated sites", Sys.time(), nrow(hyper_kept), nrow(hypo_kept)))
assertthat::assert_that(nrow(info_sites) > 0, msg = "PAMES did not find any informative site.")

purity_df <- get_purity(cbind(tumor_table, control_table), info_sites)
purity_df$Condition <- ifelse(purity_df$Sample %in% mets, "metastasis", "Normal")
dir.create(dirname(out_purity), recursive = TRUE, showWarnings = FALSE)
write.table(purity_df, out_purity, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(info_sites, out_sites, sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("[%s] Done. Purity -> %s (%d sites -> %s)", Sys.time(), out_purity, nrow(info_sites), out_sites))
