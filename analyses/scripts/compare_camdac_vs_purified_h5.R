# =============================================================================
# compare_camdac_vs_purified_h5.R
# Per CPGEA Tumour: CAMDAC package output (pure.csv.gz, m_t) vs the purified
# beta in integrated_wgbs_purified_pames.h5 (purify_tumour_NT.R, same CAMDAC
# equation applied to the integrated betas), on a fixed random sample of
# integrated CpG rows (blocks of consecutive rows).
#   - share of sampled CpGs that CAMDAC kept (the rest were imputed in the
#     integration or failed CAMDAC's effective-coverage filter)
#   - among kept CpGs with identical bulk inputs: agreement of purified betas
#   - over all kept CpGs: correlation / median |diff|
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(GenomicRanges); library(rhdf5); library(parallel) })

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
h5_dir   <- args$h5_dir
camdac   <- args$camdac_dir
n_blocks <- as.integer(if (is.null(args$n_blocks)) 100 else args$n_blocks)
n_cores  <- as.integer(if (is.null(args$n_cores)) 8 else args$n_cores)
out_tsv  <- args$out_tsv

set.seed(1)
rr <- readRDS(file.path(h5_dir, "integrated_wgbs_rowRanges.rds"))
cd <- as.data.frame(readRDS(file.path(h5_dir, "integrated_wgbs_colData.rds")))
tumours <- rownames(cd)[cd$condition == "Tumour"]
normals <- rownames(cd)[cd$condition == "Normal"]

st   <- sort(sample(seq(1, length(rr) - 999, by = 1000), n_blocks))
rows <- as.vector(outer(0:999, st, "+"))
key  <- data.table(chrom = sub("chr", "", as.character(seqnames(rr[rows]))), start = start(rr[rows]))

read_blocks <- function(f, col) unlist(lapply(st, function(s)
    h5read(file.path(h5_dir, f), "beta", index = list(s:(s + 999), col))))

one <- function(t) {
    j  <- match(t, rownames(cd))
    jn <- match(normals[match(cd[t, "patient"], cd[normals, "patient"])], rownames(cd))
    k  <- copy(key)[, `:=`(mine = read_blocks("integrated_wgbs_purified_pames.h5", j),
                          bulk = read_blocks("integrated_wgbs.h5", j),
                          bn   = read_blocks("integrated_wgbs.h5", jn))]
    p <- fread(file.path(camdac, t, "Methylation", t, sprintf("%s.%s.pure.csv.gz", t, t)),
               select = c("chrom", "start", "m", "m_n", "m_t"))
    p[, chrom := as.character(chrom)]
    x <- merge(k, p, by = c("chrom", "start"))
    same <- abs(x$bulk - x$m) < 1e-5 & abs(x$bn - x$m_n) < 1e-5
    d <- x$mine - x$m_t
    data.table(sample = t, n_sampled = nrow(k), frac_kept_by_camdac = nrow(x) / nrow(k),
               frac_same_inputs = mean(same),
               max_absdiff_same_inputs = if (any(same)) max(abs(d[same])) else NA_real_,
               cor_all_kept = cor(x$mine, x$m_t), median_absdiff_all_kept = median(abs(d)),
               frac_identical_all_kept = mean(abs(d) < 1e-5))
}

res <- rbindlist(mclapply(tumours, function(t) tryCatch(one(t), error = function(e)
    data.table(sample = t, error = conditionMessage(e))), mc.cores = n_cores), fill = TRUE)
fwrite(res, out_tsv, sep = "\t")

if ("error" %in% names(res)) print(res[!is.na(error), .(sample, error)])
ok <- res[is.na(if ("error" %in% names(res)) error else NA)]
cat(sprintf("%d / %d Tumours compared | %s CpGs sampled each\n", nrow(ok), length(tumours),
            format(ok$n_sampled[1], big.mark = ",")))
cat(sprintf("kept by CAMDAC: median %.1f%% (range %.1f-%.1f%%)\n",
            100 * median(ok$frac_kept_by_camdac), 100 * min(ok$frac_kept_by_camdac), 100 * max(ok$frac_kept_by_camdac)))
cat(sprintf("identical bulk inputs among kept: median %.1f%% (min %.1f%%)\n",
            100 * median(ok$frac_same_inputs), 100 * min(ok$frac_same_inputs)))
cat(sprintf("same inputs -> max |purified diff| across all Tumours: %.1e\n", max(ok$max_absdiff_same_inputs, na.rm = TRUE)))
cat(sprintf("all kept CpGs: identical (<1e-5) median %.2f%% (min %.2f%%) | correlation min %.5f\n",
            100 * median(ok$frac_identical_all_kept), 100 * min(ok$frac_identical_all_kept), min(ok$cor_all_kept)))
