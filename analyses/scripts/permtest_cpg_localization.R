#!/usr/bin/env Rscript
# =============================================================================
# permtest_cpg_localization.R
# Region-level permutation test (regioneR) of the genomic localization of the
# selected CpGs. Inputs from test_cpg_localization.py:
#   region_points.tsv.gz         one point per distinct selected region
#                                (its middle CpG; clusters count once)
#   candidate_categories.tsv.gz  universe: candidate CpGs of the kept models,
#                                with their exclusive category
# Random sets: as many candidate CpGs as regions, drawn per chromosome
# (resampleRegions); statistic: points in each category. Reports observed,
# permuted mean, fold, z-score and empirical p per category.
#
# Usage: Rscript permtest_cpg_localization.R <in_dir> <ntimes> <cores> [points] [out]
#   points: default region_points.tsv.gz; selected_cpg_points.tsv.gz = every selected CpG
#   out:    default regioner_permtest.tsv
# =============================================================================
suppressMessages({library(regioneR); library(data.table); library(GenomicRanges)})
args <- commandArgs(trailingOnly = TRUE)
in_dir <- args[1]; ntimes <- as.integer(args[2]); cores <- as.integer(args[3])
points_file <- if (length(args) >= 4) args[4] else "region_points.tsv.gz"
out_file <- if (length(args) >= 5) args[5] else "regioner_permtest.tsv"
set.seed(1)

u <- fread(file.path(in_dir, "candidate_categories.tsv.gz"))
p <- fread(file.path(in_dir, points_file))
cats <- fread(file.path(in_dir, "categories.tsv"))
gr <- function(d) GRanges(d$chr, IRanges(d$pos, width = 1))
U <- gr(u); A <- gr(p)
cat("universe", length(U), "CpGs; points", length(A), "\n")

B <- lapply(cats$cat, function(k) gr(u[cat == k]))
ev <- lapply(seq_along(B), function(i) { b <- B[[i]]; function(A, ...) numOverlaps(A, b, count.once = TRUE) })
names(ev) <- cats$category

pt <- permTest(A = A, ntimes = ntimes, randomize.function = resampleRegions, universe = U,
               per.chromosome = TRUE, evaluate.function = ev,
               force.parallel = TRUE, mc.cores = cores, verbose = FALSE)

res <- rbindlist(lapply(names(pt), function(k) {
  x <- pt[[k]]
  data.table(category = k, n_points = length(A), observed = x$observed,
             permuted_mean = mean(x$permuted), permuted_sd = sd(x$permuted),
             fold = x$observed / mean(x$permuted), zscore = x$zscore,
             alternative = x$alternative, p_empirical = x$pval, ntimes = x$ntimes)
}))
res[, q_BH := p.adjust(p_empirical, "BH")]
fwrite(res, file.path(in_dir, out_file), sep = "\t")
print(res)
