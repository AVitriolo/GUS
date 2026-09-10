#!/usr/bin/env Rscript
# ============================================================
# go_enrichment_methyldriver.R
# Systematic GO (BP) enrichment of significant hypo- and hyper-methylated
# genes from every MethylDriver result, with dotplots.
#
# For each *.rds in --results_dir:
#   sig = full_results with q_value < q
#   hypo  = genes with home_region_methyl_change < 0
#   hyper = genes with home_region_methyl_change > 0
#   enrichGO(BP) on each, dotplot, and a combined results table.
#
# Usage:
#   Rscript go_enrichment_methyldriver.R \
#     --results_dir=results/methyldriver \
#     --pattern='_hg38_v38_pc_50000_25_10_3.rds' \
#     --out_dir=results/go \
#     --q=0.05
# ============================================================
suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(ggplot2); library(dplyr)
})

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
results_dir <- args$results_dir
pattern     <- if (!is.null(args$pattern)) args$pattern else "\\.rds$"
out_dir     <- args$out_dir
q_cut       <- if (!is.null(args$q)) as.numeric(args$q) else 0.05
ont         <- if (!is.null(args$ont)) args$ont else "BP"
show_n      <- if (!is.null(args$show_n)) as.integer(args$show_n) else 20

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "plots"), recursive = TRUE, showWarnings = FALSE)

rds_files <- list.files(results_dir, pattern = pattern, full.names = TRUE)
# keep only the per-comparison result objects (exclude full_*.tsv etc.)
rds_files <- rds_files[grepl("\\.rds$", rds_files) &
                       !grepl("neighbors|featmatrix", basename(rds_files))]
cat(sprintf("MethylDriver result files: %d\n", length(rds_files)))

# comparison label from filename: strip the trailing _hg38_..._CV.rds
comp_of <- function(f) sub("_hg38.*$", "", basename(f))

run_go <- function(genes, label) {
    genes <- unique(genes[!is.na(genes) & nzchar(genes)])
    if (length(genes) < 5) {
        cat(sprintf("  [%s] only %d genes - skipping\n", label, length(genes)))
        return(NULL)
    }
    ego <- tryCatch(
        enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                 ont = ont, pAdjustMethod = "BH",
                 pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE),
        error = function(e) { cat(sprintf("  [%s] enrichGO error: %s\n", label, conditionMessage(e))); NULL })
    if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
        cat(sprintf("  [%s] no enriched terms\n", label)); return(NULL)
    }
    cat(sprintf("  [%s] %d genes -> %d enriched GO:%s terms\n",
                label, length(genes), nrow(as.data.frame(ego)), ont))

    # dotplot
    p <- dotplot(ego, showCategory = show_n) +
        ggtitle(label) + theme(plot.title = element_text(hjust = 0.5))
    ggsave(file.path(out_dir, "plots", paste0("GO_", ont, "_", label, ".pdf")),
           p, width = 8, height = 8)

    # table
    df <- as.data.frame(ego)
    df$set <- label
    df
}

all_tables <- list()
for (f in rds_files) {
    comp <- comp_of(f)
    cat(sprintf("== %s ==\n", comp))
    res <- readRDS(f)
    fr  <- res$full_results
    sig <- fr[!is.na(fr$q_value) & fr$q_value < q_cut, ]

    all_genes <- unique(sig$gene[!is.na(sig$gene)])   # hypo + hyper together
    t_all <- run_go(all_genes, comp)
    all_tables <- c(all_tables, list(t_all))
}


# combined table across all comparisons/directions
all_tables <- all_tables[!sapply(all_tables, is.null)]
if (length(all_tables) > 0) {
    combined <- dplyr::bind_rows(all_tables)
    out_tsv <- file.path(out_dir, sprintf("GO_%s_all_comparisons.tsv", ont))
    write.table(combined, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
    cat(sprintf("\ncombined table -> %s (%d rows)\n", out_tsv, nrow(combined)))
} else {
    cat("\nno enriched terms in any set\n")
}