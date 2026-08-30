# Concatenate all selected features (cluster names) across transcripts
# and emit a BED of their genomic coordinates.
#
# A selected feature is "<TxID>_<cluster_id>". Its coordinates come from
# the transcript's cluster map (CpGID = chr_pos  ->  cluster_id): the
# metaCpG spans from the min to the max CpG position in that cluster.
#
# Output BED columns: chr  start  end  name  n_cpgs
#   name = <TxID>_<cluster_id>
#   span mode (default): one interval per selected cluster (min..max CpG)
#   percpg mode: one line per constituent CpG of every selected cluster


args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
sel_feats_dir <- args$sel_feats_dir
map_dir       <- args$map_dir
out_bed       <- args$out_bed
mode          <- if (!is.null(args$mode)) args$mode else "span"
stopifnot(mode %in% c("span", "percpg"))

sel_files <- list.files(sel_feats_dir, full.names = TRUE)
sel_files <- sel_files[!file.info(sel_files)$isdir]
cat(sprintf("sel_feats files: %d\n", length(sel_files)))

# cache cluster maps by TxID (each map read at most once)
map_cache <- new.env(parent = emptyenv())

get_map <- function(txid) {
    if (!is.null(map_cache[[txid]])) return(map_cache[[txid]])
    # find the map file for this TxID (path has other wildcards, ends in _<TxID>_clusters_map.txt)
    hit <- list.files(map_dir,
                      pattern = sprintf("_%s_clusters_map\\.txt$", txid),
                      full.names = TRUE)
    if (length(hit) == 0) { map_cache[[txid]] <- NA; return(NA) }
    m <- tryCatch(read.table(hit[1], header = FALSE, sep = "\t",
                             stringsAsFactors = FALSE,
                             col.names = c("CpGID", "cluster_id")),
                  error = function(e) NULL)
    if (is.null(m) || nrow(m) == 0) { map_cache[[txid]] <- NA; return(NA) }
    # parse chr and pos from CpGID = "chr_pos"
    m$chr <- sub("_[0-9]+$", "", m$CpGID)
    m$pos <- as.integer(sub(".*_", "", m$CpGID))
    map_cache[[txid]] <- m
    m
}

rows <- list()
n_feat <- 0L
n_missing_map <- 0L
n_missing_cluster <- 0L

for (f in sel_files) {
    feats <- readLines(f, warn = FALSE)
    feats <- feats[nzchar(feats)]
    if (length(feats) == 0) next

    for (feat in feats) {
        n_feat <- n_feat + 1L
        # feat = "<TxID>_<cluster_id>"  ->  split on the LAST underscore
        txid <- sub("_[^_]+$", "", feat)
        cid  <- sub(".*_", "", feat)

        m <- get_map(txid)
        if (length(m) == 1 && is.na(m)) { n_missing_map <- n_missing_map + 1L; next }

        sel <- m[as.character(m$cluster_id) == cid, , drop = FALSE]
        if (nrow(sel) == 0) { n_missing_cluster <- n_missing_cluster + 1L; next }

        if (mode == "span") {
            rows[[length(rows) + 1L]] <- data.frame(
                chr    = sel$chr[1],
                start  = min(sel$pos),
                end    = max(sel$pos) + 1L,      # BED end is exclusive
                name   = feat,
                n_cpgs = nrow(sel),
                stringsAsFactors = FALSE
            )
        } else {  # percpg
            rows[[length(rows) + 1L]] <- data.frame(
                chr       = sel$chr,
                start     = sel$pos,
                end       = sel$pos + 1L,
                CpGID     = sel$CpGID,      # chr_pos, matches the H5 CpGID
                clusterID = feat,           # <TxID>_<cluster_id>
                stringsAsFactors = FALSE
            )
        }
    }
}
if (length(rows) == 0) stop("No selected features resolved to coordinates.")

bed <- do.call(rbind, rows)

# sort by chr then start
chr_rank <- function(c) {
    n <- sub("^chr", "", c)
    suppressWarnings(num <- as.integer(n))
    ifelse(is.na(num), 100L + match(n, c("X", "Y", "M", "MT")), num)
}
bed <- bed[order(chr_rank(bed$chr), bed$start), ]

# dedup on whichever columns this mode actually has
dedup_cols <- if (mode == "percpg") {
    c("chr", "start", "end", "CpGID", "clusterID")
} else {
    c("chr", "start", "end", "name")
}
bed <- bed[!duplicated(bed[, dedup_cols]), ]

dir.create(dirname(out_bed), recursive = TRUE, showWarnings = FALSE)
write.table(bed, out_bed, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)