#!/usr/bin/env Rscript
# =============================================================================
# sel_feats_to_regions_w20.R
# Regions for methyldriver / gimme from the per-CpG (no clustering) models.
#   - selected CpGs of the SAME transcript that are <= --gap bp apart are
#     chained into one region (A-B <= gap and B-C <= gap -> A,B,C one region)
#   - a region spans first CpG - flank .. last CpG + flank (BED, 0-based
#     half-open: [pos - flank, pos + flank) for a lone CpG = 2 * flank bp)
#   - models kept: cv_r2 > mean cv_r2 of all models (--keep=all: every model with features)
#   - methyldriver unit: one region per distinct CpG set; transcripts whose
#     regions are exactly the same are collapsed into it (TxIDs listed)
#
# Outputs (in --out_dir, suffixed _<tag>):
#   features_*.tsv      TxID CpGID pos cv_r2 regionID coordID  (all models)
#   regions_*.tsv       one row per transcript region            (all models)
#   coords_all_*.bed    distinct windows of all regions (bigWig averages)
#   md_regions_*.tsv    methyldriver regions (models above the mean)
#   md_percpg_*.bed     chr start end CpGID mdID (cluster_delta_methylation.R)
#   gimme_coords_*.bed  distinct windows of the methyldriver regions
#   cvr2_*.tsv          TxID cv_r2 kept
# =============================================================================
suppressPackageStartupMessages(library(data.table))
options(scipen = 999)

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
sel_dir  <- args$sel_feats_dir
perf_dir <- args$perf_dir
out_dir  <- args$out_dir
tag      <- args$tag
gap      <- as.integer(if (is.null(args$gap)) 20 else args$gap)
flank    <- as.integer(if (is.null(args$flank)) 10 else args$flank)
keep_all <- !is.null(args$keep) && args$keep == "all"   # --keep=all: every model (e.g. features already filtered by stability)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out <- function(stem, ext) file.path(out_dir, sprintf("%s_%s.%s", stem, tag, ext))

# ---- selected features: one line per feature, TxID_chr_pos ----
feat <- fread(cmd = sprintf("find %s -type f -exec cat {} +", shQuote(sel_dir)), header = FALSE, col.names = "f")
feat <- feat[nzchar(f)]
feat[, `:=`(TxID = sub("_.*$", "", f), CpGID = sub("^[^_]*_", "", f))]
feat[, `:=`(chr = sub("_[0-9]+$", "", CpGID), pos = as.integer(sub(".*_", "", CpGID)), f = NULL)]
feat <- unique(feat)
stopifnot(feat[, uniqueN(chr), by = TxID][, all(V1 == 1)])
cat(sprintf("selected features: %d | transcripts: %d | distinct CpGs: %d\n",
            nrow(feat), uniqueN(feat$TxID), uniqueN(feat$CpGID)))

# ---- cv_r2 per model (file name ends with _<TxID>) ----
perf <- fread(cmd = sprintf("find %s -type f -exec awk -v OFS='\\t' '$1==\"cv_r2\" {n=split(FILENAME,a,\"_\"); print a[n], $2}' {} +",
                            shQuote(perf_dir)), header = FALSE, col.names = c("TxID", "cv_r2"))
stopifnot(!anyDuplicated(perf$TxID), all(feat$TxID %in% perf$TxID))
mean_cv <- mean(perf$cv_r2)
perf[, kept := if (keep_all) TxID %in% feat$TxID else cv_r2 > mean_cv]
cat(sprintf("models: %d | mean cv_r2 %.4f | kept (%s): %d\n", nrow(perf), mean_cv,
            if (keep_all) "--keep=all: every model with features" else "cv_r2 > mean", sum(perf$kept)))
fwrite(perf, out("cvr2", "tsv"), sep = "\t")

# ---- chain CpGs of the same transcript ----
setorder(feat, TxID, pos)
feat[, grp := cumsum(c(TRUE, diff(pos) > gap)), by = TxID]
feat[, `:=`(start = min(pos) - flank, end = max(pos) + flank), by = .(TxID, grp)]
feat[, `:=`(regionID = paste(TxID, chr, start, end, sep = "_"), coordID = paste(chr, start, end, sep = "_"))]
feat[perf, cv_r2 := i.cv_r2, on = "TxID"]
fwrite(feat[, .(TxID, CpGID, pos, cv_r2, regionID, coordID)], out("features", "tsv"), sep = "\t")

reg <- feat[, .(TxID = TxID[1], chr = chr[1], start = start[1], end = end[1], n_CpGs = .N,
                CpGs = paste(CpGID, collapse = ","), cv_r2 = cv_r2[1], coordID = coordID[1]), by = regionID]
cat(sprintf("transcript regions: %d (%d with >1 CpG) | width median %d, max %d bp\n",
            nrow(reg), sum(reg$n_CpGs > 1), as.integer(median(reg$end - reg$start)), max(reg$end - reg$start)))
fwrite(reg, out("regions", "tsv"), sep = "\t")

bed_of <- function(d) setorder(unique(d[, .(chr, start, end, coordID)]), chr, start)
fwrite(bed_of(reg), out("coords_all", "bed"), sep = "\t", col.names = FALSE)

# ---- methyldriver regions: models above the mean, one per distinct CpG set ----
kept <- reg[TxID %in% perf[kept == TRUE, TxID]]
md <- kept[, .(chr = chr[1], start = start[1], end = end[1], n_CpGs = n_CpGs[1],
               TxIDs = paste(sort(unique(TxID)), collapse = ","), n_tx = uniqueN(TxID)),
           by = .(coordID, CpGs)]
# same window, different CpG sets (a transcript skipped an inner CpG): kept apart
md[, mdID := if (.N == 1) coordID else paste0(coordID, "_v", seq_len(.N)), by = coordID]
setcolorder(md, c("mdID", "chr", "start", "end", "n_CpGs", "n_tx", "TxIDs", "coordID", "CpGs"))
setorder(md, chr, start, end)
cat(sprintf("methyldriver regions: %d from %d transcript regions (%d shared by >1 transcript; %d windows with >1 CpG set)\n",
            nrow(md), nrow(kept), sum(md$n_tx > 1), md[, .N, by = coordID][N > 1, .N]))
fwrite(md, out("md_regions", "tsv"), sep = "\t")

percpg <- md[, .(CpGID = strsplit(CpGs, ",", fixed = TRUE)[[1]]), by = .(mdID, chr)]
percpg[, pos := as.integer(sub(".*_", "", CpGID))]
fwrite(percpg[, .(chr, start = pos - flank, end = pos + flank, CpGID, mdID)],
       out("md_percpg", "bed"), sep = "\t", col.names = FALSE)

fwrite(bed_of(md), out("gimme_coords", "bed"), sep = "\t", col.names = FALSE)
cat("done\n")
