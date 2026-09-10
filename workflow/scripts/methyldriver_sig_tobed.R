#!/usr/bin/env Rscript
# ============================================================
# methyldriver_sig_to_bed.R
# Build a per-CpG BED of MethylDriver-significant clusters.
#   1. collect significant clusters (q_value < q) from MethylDriver full tables
#   2. expand each cluster to its constituent CpGs via the cluster maps
#   3. write chr start end CpGID clusterID   (one row per CpG)
#
# Significant clusters are pooled across all --full_dir tables matching
# --pattern (union of clusters significant in ANY comparison). Use
# --full_file for a single comparison instead.
#
# Usage:
#   Rscript methyldriver_sig_to_bed.R \
#     --full_dir=results/methyldriver \
#     --pattern='^full_.*\.tsv$' \
#     --map_dir=data/clustering/clusters_maps \
#     --q=0.05 \
#     --out_bed=data/to_methyldriver/methyldriver_sig_cpgs.bed
# ============================================================
suppressPackageStartupMessages({ library(dplyr) })

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
full_dir  <- args$full_dir
full_file <- args$full_file
pattern   <- if (!is.null(args$pattern)) args$pattern else "^full_.*\\.tsv$"
map_dir   <- args$map_dir
q_cut     <- if (!is.null(args$q)) as.numeric(args$q) else 0.05
out_bed   <- args$out_bed

dir.create(dirname(out_bed), recursive = TRUE, showWarnings = FALSE)

# ---- 1. gather significant cluster IDs ----
if (!is.null(full_file)) {
    full_files <- full_file
} else {
    full_files <- list.files(full_dir, pattern = pattern, full.names = TRUE)
}
cat(sprintf("MethylDriver full tables: %d\n", length(full_files)))
stopifnot(length(full_files) > 0)

sig_clusters <- unique(unlist(lapply(full_files, function(f) {
    d <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    # the cluster id column is 'name' in MethylDriver full_results
    idcol <- if ("name" %in% names(d)) "name" else names(d)[1]
    d[[idcol]][!is.na(d$q_value) & d$q_value < q_cut]
})))
cat(sprintf("significant clusters (q < %.3g), pooled: %d\n", q_cut, length(sig_clusters)))
stopifnot(length(sig_clusters) > 0)

# ---- 2. expand to CpGs via cluster maps ----
# cluster id = "<TxID>_<cluster_id>"; map file per TxID has CpGID<TAB>cluster_id
sig_by_tx <- split(sig_clusters, sub("_[^_]+$", "", sig_clusters))  # group by TxID

map_cache <- new.env(parent = emptyenv())
get_map <- function(txid) {
    if (!is.null(map_cache[[txid]])) return(map_cache[[txid]])
    hit <- list.files(map_dir, pattern = sprintf("_%s_clusters_map\\.txt$", txid),
                      full.names = TRUE)
    m <- if (length(hit) == 0) NULL else
        tryCatch(read.table(hit[1], header = FALSE, sep = "\t",
                            stringsAsFactors = FALSE,
                            col.names = c("CpGID", "cluster_id")),
                 error = function(e) NULL)
    map_cache[[txid]] <- m
    m
}

rows <- list()
n_missing <- 0L
for (txid in names(sig_by_tx)) {
    m <- get_map(txid)
    if (is.null(m)) { n_missing <- n_missing + length(sig_by_tx[[txid]]); next }
    for (feat in sig_by_tx[[txid]]) {
        cid <- sub(".*_", "", feat)                    # trailing cluster number
        sel <- m[as.character(m$cluster_id) == cid, , drop = FALSE]
        if (nrow(sel) == 0) { n_missing <- n_missing + 1L; next }
        chr <- sub("_[0-9]+$", "", sel$CpGID)
        pos <- as.integer(sub(".*_", "", sel$CpGID))
        rows[[length(rows) + 1L]] <- data.frame(
            chr = chr, start = pos, end = pos + 1L,
            CpGID = sel$CpGID, clusterID = feat, stringsAsFactors = FALSE)
    }
}
if (n_missing > 0) cat(sprintf("clusters not resolved to CpGs: %d\n", n_missing))

bed <- do.call(rbind, rows)
bed <- bed[!duplicated(bed[, c("chr","start","end","CpGID","clusterID")]), ]

# sort chr then start
chr_rank <- function(c) { n <- sub("^chr","",c); suppressWarnings(v <- as.integer(n))
    ifelse(is.na(v), 100L + match(n, c("X","Y","M","MT")), v) }
bed <- bed[order(chr_rank(bed$chr), bed$start), ]

write.table(bed, out_bed, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
cat(sprintf("wrote %d CpG rows (%d clusters) -> %s\n",
            nrow(bed), length(unique(bed$clusterID)), out_bed))