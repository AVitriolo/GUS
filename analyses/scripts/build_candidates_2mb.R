# Per-TxID candidate CpG lists (TxID<TAB>CpGID) = all CpGs of an integrated rowRanges inside the
# TxID's +-1 Mb TSS window (tss_all_2mb.window.sorted.bed, BED: start 0-based, end inclusive).
suppressPackageStartupMessages({library(GenomicRanges); library(data.table)})
a <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
rr <- readRDS(a$rowranges); w <- fread(a$windows, header = FALSE, select = 1:4, col.names = c("chr", "start", "end", "TxID"))
if (!is.null(a$only)) w <- w[TxID %in% strsplit(a$only, ",")[[1]]]
win <- GRanges(w$chr, IRanges(w$start + 1L, w$end))
ov <- findOverlaps(win, rr)
id <- paste0(as.character(seqnames(rr)), "_", start(rr))
dt <- data.table(TxID = w$TxID[queryHits(ov)], CpGID = id[subjectHits(ov)], o = subjectHits(ov))
setorder(dt, TxID, o)
dir.create(a$out_dir, showWarnings = FALSE, recursive = TRUE)
dt[, fwrite(.SD[, .(TxID, CpGID)], file.path(a$out_dir, sprintf("hg38_v38_pc_50000_25_10_%s.bed", .BY$TxID)), sep = "\t", col.names = FALSE), by = TxID]
cat(sprintf("%d TxIDs, %d candidate lines, median %d per TxID\n", uniqueN(dt$TxID), nrow(dt), as.integer(median(dt[, .N, by = TxID]$N))))
