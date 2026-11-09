# Stability selection with the models' own selector, step 3: REGION-level stability.
# Neighbouring CpGs carry the same signal and compete for the k slots, so single-CpG selection
# frequencies are split between them. Here the gene's candidate CpGs are grouped into regions by
# genomic distance only (bumphunter::clusterMaker, consecutive CpGs <= max_gap bp apart; expression
# is not used), and a region counts as selected in a half-sample if any of its CpGs was selected.
#
# Per gene (region_genes.tsv):
#   n_regions, q_regions   number of candidate regions (p) and mean selected regions per half-sample (q)
#   phi_region             Nogueira stability index of the region selections (stabm)
#   cutoff_pfer1           region frequency cutoff giving PFER <= 1 (stabs, complementary pairs,
#                          unimodality; q rounded up), NA if not attainable
#   pct_saved_regions_stable  % of regions holding >= 1 of the model's saved selected CpGs that are stable
#   pct_saved_cpgs_stable_region  % of the model's saved selected CpGs lying in a stable region
# Per region (region_freq.tsv): coordinates, n CpGs, n saved selected CpGs, frequency, stable flag.
#
# Optional region_map = name of precomputed regions (e.g. corr0.5 -> region_map_corr0.5_shard*.tsv from
# stability_regions_corr.py: distance AND adjacent-CpG correlation); replaces the distance-only clusterMaker
# regions, outputs region_genes_<name>.tsv / region_freq_<name>.tsv. bumphunter is only needed without a region_map.
# The map may carry an 8th column chr (needed when the CpG IDs are not chr_pos, e.g. 450K probe IDs).
#
# Usage: Rscript stability_selection_regions.R <in_dir> <cand_dir> [max_gap=300] [n_pairs=100] [pfer=1] [cores=12] [region_map]

suppressPackageStartupMessages({ library(stabm); library(stabs); library(data.table); library(parallel) })

args     <- commandArgs(trailingOnly = TRUE)
in_dir   <- args[1]
cand_dir <- args[2]
max_gap  <- if (length(args) >= 3) as.integer(args[3]) else 300L
n_pairs  <- if (length(args) >= 4) as.integer(args[4]) else 100L
pfer     <- if (length(args) >= 5) as.numeric(args[5]) else 1
cores    <- if (length(args) >= 6) as.integer(args[6]) else 12L
map_name <- if (length(args) >= 7) args[7] else NA_character_
out_tag  <- if (is.na(map_name)) sprintf("gap%d", max_gap) else map_name
prefix   <- "hg38_v38_pc_50000_25_10"
B <- 2L * n_pairs

rmap <- NULL
if (!is.na(map_name)) {
    rmap <- rbindlist(lapply(list.files(in_dir, sprintf("^region_map_%s_shard.*\\.tsv$", map_name), full.names = TRUE),
                             function(f) {
                                 x <- fread(f, sep = "\t", header = FALSE)
                                 # optional 8th column chr (CpG IDs that do not carry it, e.g. 450K probe IDs)
                                 setnames(x, c("TxID", "n_regions", "CpGID", "region", "start", "end", "n_cpg", "chr")[seq_len(ncol(x))])
                             }), fill = TRUE)
    rmap <- split(rmap, by = "TxID")
}

sel  <- rbindlist(lapply(list.files(in_dir, "^subsample_selections_shard.*\\.tsv$", full.names = TRUE),
                         fread, sep = "\t", header = FALSE, col.names = c("TxID", "b", "idx"), colClasses = "character"))
cand <- rbindlist(lapply(list.files(in_dir, "^candidates_shard.*\\.tsv$", full.names = TRUE),
                         fread, sep = "\t", header = FALSE,
                         col.names = c("TxID", "n_cand", "k", "idx", "CpGID", "saved_selected")))
sel  <- split(sel$idx, sel$TxID)
cand <- split(cand, by = "TxID")

one_gene <- function(TxID) {
    sets <- sel[[TxID]]
    cg   <- cand[[TxID]]
    if (length(sets) != B || is.null(cg)) return(NULL)

    if (is.null(rmap)) {
        bed <- fread(file.path(cand_dir, sprintf("%s_%s.bed", prefix, TxID)), header = FALSE, select = 2L,
                     col.names = "CpGID")
        bed <- unique(bed)
        bed[, c("chr", "pos") := tstrsplit(CpGID, "_", fixed = TRUE)][, pos := as.integer(pos)]
        setorder(bed, pos)
        bed[, region := bumphunter::clusterMaker(chr, pos, maxGap = max_gap)]
        n_reg  <- max(bed$region)
        coords <- bed[, .(start = min(pos), end = max(pos), n_cpg = .N), by = .(chr, region)]
    } else {
        bed <- rmap[[TxID]]
        if (is.null(bed)) return(NULL)
        n_reg  <- bed$n_regions[1]
        if (!"chr" %in% names(bed) || all(is.na(bed$chr))) bed[, chr := sub("_.*", "", CpGID)]
        coords <- unique(bed[, .(chr, region, start, end, n_cpg)])
    }
    region_of <- setNames(bed$region, bed$CpGID)
    cpg_of    <- setNames(cg$CpGID, cg$idx)

    reg_sets <- lapply(strsplit(sets, ","), function(s) unique(region_of[cpg_of[s]]))
    q        <- mean(lengths(reg_sets))
    cutoff   <- tryCatch(stabsel_parameters(p = n_reg, q = ceiling(q), PFER = pfer, B = n_pairs,
                                            sampling.type = "SS", assumption = "unimodal")$cutoff,
                         error = function(e) NA_real_)

    freq <- table(unlist(reg_sets)) / B
    saved_reg <- region_of[cg$CpGID[cg$saved_selected == 1]]
    reg <- coords[region %in% c(as.integer(names(freq)), saved_reg)]
    reg[, freq := as.numeric(freq[as.character(region)])][is.na(freq), freq := 0]
    reg[, n_saved := as.integer(table(factor(saved_reg, levels = region)))]
    reg[, stable := !is.na(cutoff) & freq >= cutoff]
    reg[, TxID := TxID]

    gene <- data.table(TxID = TxID, k = cg$k[1], n_regions = n_reg, q_regions = q,
                       phi_region = stabilityNogueira(reg_sets, p = n_reg), cutoff_pfer1 = cutoff,
                       n_stable_regions = sum(reg$stable),
                       n_saved_regions = sum(reg$n_saved > 0),
                       pct_saved_regions_stable = 100 * mean(reg$stable[reg$n_saved > 0]),
                       pct_saved_cpgs_stable_region = 100 * sum(reg$n_saved[reg$stable]) / sum(reg$n_saved))
    list(gene = gene, reg = reg)
}

res <- mclapply(names(sel), one_gene, mc.cores = cores)
res <- res[!vapply(res, is.null, logical(1))]
genes <- rbindlist(lapply(res, `[[`, "gene"))
regs  <- rbindlist(lapply(res, `[[`, "reg"))
setcolorder(regs, c("TxID", "chr", "region", "start", "end", "n_cpg", "n_saved", "freq", "stable"))

fwrite(genes, file.path(in_dir, sprintf("region_genes_%s.tsv", out_tag)), sep = "\t")
fwrite(regs[order(TxID, -freq)], file.path(in_dir, sprintf("region_freq_%s.tsv", out_tag)), sep = "\t")
cat(sprintf("%d genes, regions %s; cutoff NA for %d genes\n", nrow(genes), out_tag, sum(is.na(genes$cutoff_pfer1))))
cat(sprintf("saved selected CpGs in a stable region: %.1f%% (%d / %d)\n",
            100 * sum(regs$n_saved[regs$stable]) / sum(regs$n_saved), sum(regs$n_saved[regs$stable]), sum(regs$n_saved)))
cat(sprintf("regions with >=1 saved CpG that are stable: %.1f%%\n", 100 * mean(regs$stable[regs$n_saved > 0])))
print(summary(genes[, .(n_regions, q_regions, phi_region, cutoff_pfer1, pct_saved_cpgs_stable_region)]))
