#!/usr/bin/env Rscript
# =============================================================================
# methyldriver_w20.R
# MethylDriver on the w20 regions (sel_feats_to_regions_w20.R).
#   Neighbours, once for all comparisons (helpers/neighbors_min_dist.R): PCA
#   of the chromatin feature matrix (centred, scaled, first --n_pcs PCs) + kNN,
#   as build_neighbors.R / most_similar_n_regions.R, but a candidate on the
#   same chromosome closer than --min_dist bp (gap between the two windows;
#   overlapping and identical windows included) is skipped and the next most
#   similar region is taken, so a region is not compared with its own locus.
#   Reused if neighbors_<tag>.rds exists (e.g. built on the HPC by
#   build_neighbors_w20.R).
#   Then per comparison: z of the region's delta against its --n neighbours'
#   deltas (helpers/MethylDriver.R, unchanged).
# A region collapses the transcripts whose regions are exactly the same, so
# its genes are all the genes of those transcripts (comma-separated).
# =============================================================================
set.seed(123); options(scipen = 999)
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel) })
source("workflow/scripts/helpers/MethylDriver.R")
source("analyses/scripts/helpers/neighbors_min_dist.R")

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
fm_path    <- args$featmatrix
md_path    <- args$md_regions
deltas     <- strsplit(strsplit(args$deltas, ",")[[1]], "=")       # comparison=path,...
deltas     <- setNames(sapply(deltas, `[`, 2), sapply(deltas, `[`, 1))
conv_path  <- args$conversion_table
out_dir    <- args$out_dir
tag        <- args$tag
n_nb       <- as.integer(if (is.null(args$n)) 100 else args$n)
n_pcs      <- as.integer(if (is.null(args$n_pcs)) 20 else args$n_pcs)
min_dist   <- as.integer(if (is.null(args$min_dist)) 1000 else args$min_dist)
q_cutoff   <- as.numeric(if (is.null(args$q)) 0.05 else args$q)
num_cores  <- as.integer(if (is.null(args$num_cores)) 8 else args$num_cores)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out <- function(stem, ext) file.path(out_dir, sprintf("%s_%s.%s", stem, tag, ext))

md <- fread(md_path)
N  <- nrow(md)

# ---- neighbours ----
nb_file <- out("neighbors", "rds")
if (file.exists(nb_file)) {
    closest_regions_list <- readRDS(nb_file)
    stopifnot(identical(names(closest_regions_list), md$mdID))
    message("neighbours: reusing ", nb_file)
} else {
    fm <- readRDS(fm_path)
    x  <- fm[md$coordID, , drop = FALSE]; rm(fm); invisible(gc())
    rownames(x) <- md$mdID
    message(sprintf("feature matrix: %d regions x %d features", nrow(x), ncol(x)))
    closest_regions_list <- most_similar_n_regions_min_dist(x, md$chr, md$start, md$end, n = n_nb, n_pcs = n_pcs,
                                                            min_dist = min_dist, num_cores = num_cores)
    rm(x); invisible(gc())
    saveRDS(closest_regions_list, nb_file)
}

# ---- genes of each region ----
conv <- read.table(conv_path, header = TRUE, stringsAsFactors = FALSE)
long <- md[, .(TxID = strsplit(TxIDs, ",", fixed = TRUE)[[1]]), by = mdID]
long[, gene := conv$gene_name[match(TxID, conv$transcript_id)]]
genes <- long[, .(genes = paste(unique(na.omit(gene)), collapse = ",")), by = mdID][md$mdID, on = "mdID"]$genes
split_genes <- function(g) unique(unlist(strsplit(g, ",", fixed = TRUE)))

# ---- per comparison ----
for (cmp in names(deltas)) {
    d <- fread(deltas[[cmp]])
    delta <- setNames(d$delta, d$clusterID)[md$mdID]
    if (anyNA(delta)) stop(sprintf("%s: %d regions without a delta", cmp, sum(is.na(delta))))

    res <- MethylDriver(home_region_methyl_change = delta, regions_genes = genes,
                        closest_regions_list = closest_regions_list, q_value_cutoff = q_cutoff)
    full <- as.data.table(res$full_results)
    full[, gene := as.character(gene)]
    full <- md[, .(mdID, chr, start, end, n_CpGs, n_tx, TxIDs)][full, on = c(mdID = "name")]
    setnames(full, "mdID", "name")
    setorder(full, q_value)
    res$full_results <- as.data.frame(full)
    res$hypermethylated_genes <- split_genes(full[q_value < q_cutoff & home_region_methyl_change > 0, gene])
    res$hypomethylated_genes  <- split_genes(full[q_value < q_cutoff & home_region_methyl_change < 0, gene])
    saveRDS(res, out(cmp, "rds"))
    fwrite(full, out(paste0("full_", cmp), "tsv"), sep = "\t")

    full[, significance := fifelse(q_value < q_cutoff & home_region_methyl_change > 0, "Hyper",
                           fifelse(q_value < q_cutoff & home_region_methyl_change < 0, "Hypo", "NS"))]
    top <- full[order(p_value)][seq_len(min(55, .N))]
    pdf(out(paste0("volcano_", cmp), "pdf"), width = 12, height = 12)
    print(ggplot(full, aes(home_region_methyl_change, -log10(p_value), color = significance)) +
      geom_point(alpha = .7) +
      ggrepel::geom_text_repel(data = top, aes(label = gene), color = "black", size = 6, max.overlaps = 100) +
      scale_color_manual(values = c(Hypo = "steelblue3", Hyper = "#F8766B", NS = "grey")) +
      labs(x = "Methylation Change", y = "-log10(p-value)", title = "MethylDriver") +
      theme_minimal())
    dev.off()
    message(sprintf("%s: %d significant regions | %d hyper, %d hypo genes", cmp,
                    sum(full$q_value < q_cutoff), length(res$hypermethylated_genes), length(res$hypomethylated_genes)))
}
