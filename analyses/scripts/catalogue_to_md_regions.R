# =============================================================================
# catalogue_to_md_regions.R
# Saturation catalogue (catalogue_cpgs.tsv.gz: cpg = HM450 probe ID or chr_pos,
# n_cohorts, first_step, first_cohort, cohorts) -> methyldriver regions, in the
# format of sel_feats_to_regions_w20.R (md_regions + md_percpg BED).
#   positions : probe -> manifest CpG_beg + 1 (1-based C of the CpG, the WGBS
#               chr_pos convention); chr_pos kept; every position must be an
#               hg38 reference CpG (checked against --se_dir rows)
#   regions   : catalogue CpGs <= --gap bp apart chained (all catalogue CpGs
#               together, no gene grouping), window first - flank .. last + flank
#   genes     : with --pairs (catalogue_pairs.tsv.gz of saturation_catalogue.py
#               --pairs: gene, cpg, ...), the genes whose models selected the
#               region's CpGs (model_genes, any cohort); TxIDs = one transcript
#               per model gene (a protein-coding TSS of --tss_bed when the gene
#               has one, else its first transcript in --conversion_table), so
#               methyldriver_w20.R labels the region with its model genes.
#               Always also: nearest protein-coding TSS (gencode v38) as
#               nearest_TxID / nearest_gene, and dist_tss (bp, 0 = overlapping).
#               Without --pairs, or for regions whose model genes are not in
#               the conversion table, TxIDs = the nearest TSS (gene_source).
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(SummarizedExperiment); library(HDF5Array)
})
args  <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
gap   <- as.integer(if (is.null(args$gap)) 20 else args$gap)
flank <- as.integer(if (is.null(args$flank)) 10 else args$flank)
dir.create(dirname(args$out_regions), recursive = TRUE, showWarnings = FALSE)

cat_ <- fread(args$catalogue)
man  <- fread(args$manifest, select = c("Probe_ID", "CpG_chrm", "CpG_beg"))
is_probe <- grepl("^(cg[0-9]|ch\\.)", cat_$cpg)
m <- man[match(cat_$cpg[is_probe], Probe_ID)]
cat(sprintf("catalogue: %d entries (%d probes, %d chr_pos); probes without manifest position: %d\n",
            nrow(cat_), sum(is_probe), sum(!is_probe), sum(is.na(m$CpG_beg))))
pos <- rbind(data.table(entry = cat_$cpg[is_probe], chr = m$CpG_chrm, pos = m$CpG_beg + 1L),
             data.table(entry = cat_$cpg[!is_probe], chr = sub("_.*", "", cat_$cpg[!is_probe]),
                        pos = as.integer(sub(".*_", "", cat_$cpg[!is_probe]))))
pos <- pos[!is.na(pos) & chr %in% paste0("chr", c(1:22, "X", "Y"))]
pos <- merge(pos, cat_[, .(entry = cpg, n_cohorts, cohorts)], by = "entry")
pos[, CpGID := paste(chr, pos, sep = "_")]

# every position must be an hg38 reference CpG (same rows as the WGBS SEs)
se <- loadHDF5SummarizedExperiment(args$se_dir)
ref <- paste(as.character(seqnames(se)), start(se), sep = "_")
ok <- pos$CpGID %in% ref
cat(sprintf("positions on an hg38 reference CpG: %d of %d (dropped %d)\n", sum(ok), nrow(pos), sum(!ok)))
pos <- pos[ok]
# a probe and a chr_pos entry can be the same CpG: one row per CpG
cpg <- pos[, .(entries = paste(entry, collapse = ","), n_cohorts = max(n_cohorts),
               cohorts = paste(unique(unlist(strsplit(cohorts, ","))), collapse = ",")), by = .(CpGID, chr, pos)]
setorder(cpg, chr, pos)
cat(sprintf("distinct catalogue CpGs: %d\n", nrow(cpg)))

cpg[, grp := cumsum(c(TRUE, diff(pos) > gap)), by = chr]
reg <- cpg[, .(start = min(pos) - flank, end = max(pos) + flank, n_CpGs = .N, CpGs = paste(CpGID, collapse = ","),
               n_cohorts = max(n_cohorts), cohorts = paste(unique(unlist(strsplit(cohorts, ","))), collapse = ",")),
           by = .(chr, grp)]
reg[, `:=`(mdID = sprintf("%s_%d_%d", chr, start, end), grp = NULL)]
reg[, coordID := mdID]

tss <- fread(args$tss_bed, header = FALSE, col.names = c("chr", "s", "e", "TxID"))
tss_gr <- GRanges(tss$chr, IRanges(tss$e, width = 1))
reg_gr <- GRanges(reg$chr, IRanges(reg$start + 1L, reg$end))
nn <- distanceToNearest(reg_gr, tss_gr, ignore.strand = TRUE)
reg[queryHits(nn), `:=`(nearest_TxID = tss$TxID[subjectHits(nn)], dist_tss = mcols(nn)$distance)]
reg[is.na(nearest_TxID), `:=`(nearest_TxID = "", dist_tss = NA_integer_)]
conv <- fread(args$conversion_table, header = FALSE, skip = 1, select = 5:6, col.names = c("TxID", "gene"))
reg[, nearest_gene := conv$gene[match(nearest_TxID, conv$TxID)]]
reg[, `:=`(TxIDs = nearest_TxID, model_genes = "", gene_source = "nearest_tss")]
if (!is.null(args$pairs)) {
  pr <- fread(args$pairs, select = c("gene", "cpg"))
  pr <- unique(merge(pr, pos[, .(cpg = entry, CpGID)], by = "cpg")[, .(CpGID, gene)])
  # gene -> one transcript: a protein-coding TSS of the gene if any, else its first transcript
  tx_of <- conv[, .(TxID = if (any(TxID %in% tss$TxID)) TxID[TxID %in% tss$TxID][1] else TxID[1]), by = gene]
  rg <- merge(reg[, .(CpGID = unlist(strsplit(CpGs, ",", fixed = TRUE))), by = mdID], pr, by = "CpGID")
  rg <- rg[, .(model_genes = paste(sort(unique(gene)), collapse = ","),
               model_TxIDs = paste(na.omit(unique(tx_of$TxID[match(gene, tx_of$gene)])), collapse = ",")), by = mdID]
  reg[rg, on = "mdID", model_genes := i.model_genes]
  reg[rg[model_TxIDs != ""], on = "mdID", `:=`(TxIDs = i.model_TxIDs, gene_source = "model")]
  cat(sprintf("model genes: %d of %d regions (%d labelled by them; %d fall back to the nearest TSS)\n",
              sum(reg$model_genes != ""), nrow(reg), sum(reg$gene_source == "model"), sum(reg$gene_source != "model")))
}
reg[, n_tx := lengths(strsplit(TxIDs, ",", fixed = TRUE))]
cat(sprintf("regions: %d (CpGs per region: median %d, max %d); distance to nearest TSS: median %s bp\n",
            nrow(reg), as.integer(median(reg$n_CpGs)), max(reg$n_CpGs), format(median(reg$dist_tss, na.rm = TRUE), big.mark = ",")))

setcolorder(reg, c("mdID", "chr", "start", "end", "n_CpGs", "n_tx", "TxIDs", "coordID", "CpGs", "model_genes", "gene_source",
                   "nearest_TxID", "nearest_gene", "dist_tss", "n_cohorts", "cohorts"))
fwrite(reg, args$out_regions, sep = "\t")
percpg <- merge(cpg[, .(CpGID, chr, pos)], reg[, .(mdID, CpGID = strsplit(CpGs, ",", fixed = TRUE)), by = mdID][, .(CpGID = unlist(CpGID)), by = mdID], by = "CpGID")
fwrite(percpg[, .(chr, start = pos - flank, end = pos + flank, CpGID, mdID)][order(chr, start)], args$out_percpg,
       sep = "\t", col.names = FALSE)
cat("written", args$out_regions, "and", args$out_percpg, "\n")
