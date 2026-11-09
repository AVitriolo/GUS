#!/usr/bin/env Rscript
# =============================================================================
# go_methyldriver_sites.R
# Metastatic-site differences from MethylDriver (methyldriver_w20.R): genes of the significant regions (q < --q) of
# the three site-vs-site comparisons, their Venn diagram, and GO BP enrichment of the site-specific genes.
#   bone_lymph / bone_liver / lymph_liver : D_Bone_vs_D_Lymph, D_Bone_vs_D_Liver, D_Lymph_vs_D_Liver
#   site_<s>       : significant in both comparisons involving site s and not in the third (s differs from both)
#   unique_<a_b>   : significant only in the a-vs-b comparison;  shared_all_three : in all three
#   all_site_genes : union of the three comparisons
# GO: enrichGO (BP, BH), universe = every gene with a tested MethylDriver region (the gene background of the run);
# sets with < --min_genes genes are not tested. Same set definitions as results/methyldriver_2mb_full_cv3r5/go_sites.
# Usage: Rscript go_methyldriver_sites.R --res_dir=DIR --tag=TAG --out_dir=DIR [--q=0.05] [--min_genes=5]
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(clusterProfiler); library(org.Hs.eg.db); library(ggplot2) })

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
res_dir <- args$res_dir; tag <- args$tag; out_dir <- args$out_dir
q_cut <- if (is.null(args$q)) 0.05 else as.numeric(args$q)
min_genes <- if (is.null(args$min_genes)) 5L else as.integer(args$min_genes)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cmp <- c(bone_lymph = "D_Bone_vs_D_Lymph", bone_liver = "D_Bone_vs_D_Liver", lymph_liver = "D_Lymph_vs_D_Liver")
split_genes <- function(x) unique(unlist(strsplit(x[!is.na(x) & nzchar(x)], ",")))
full <- lapply(cmp, function(c) fread(file.path(res_dir, sprintf("full_%s_%s.tsv", c, tag))))
universe <- split_genes(full[[1]]$gene)
sig <- lapply(full, function(d) split_genes(d[q_value < q_cut, gene]))
cat(sprintf("universe: %d genes | significant genes: %s\n", length(universe),
            paste(sprintf("%s %d", names(sig), lengths(sig)), collapse = ", ")))

# ---- Venn membership ----
genes <- sort(unique(unlist(sig)))
mem <- data.table(gene = genes, bone_lymph = genes %in% sig$bone_lymph, bone_liver = genes %in% sig$bone_liver,
                  lymph_liver = genes %in% sig$lymph_liver)
mem[, region := fcase(bone_lymph & bone_liver & lymph_liver, "shared_all_three",
                      bone_lymph & bone_liver, "site_bone", bone_liver & lymph_liver, "site_liver",
                      bone_lymph & lymph_liver, "site_lymph_node", bone_lymph, "unique_bone_lymph",
                      bone_liver, "unique_bone_liver", default = "unique_lymph_liver")]
fwrite(mem, file.path(out_dir, "venn_membership.tsv"), sep = "\t")
sets <- c(split(mem$gene, mem$region), list(all_site_genes = genes))
for (s in names(sets)) writeLines(sets[[s]], file.path(out_dir, paste0(s, "_genes.txt")))
print(mem[, .N, by = region][order(region)])

# ---- Venn diagram (3 circles) ----
circ <- function(x0, y0, r = 1.5, n = 200) data.table(x = x0 + r * cos(seq(0, 2 * pi, length.out = n)),
                                                     y = y0 + r * sin(seq(0, 2 * pi, length.out = n)))
cc <- list(bone_lymph = c(-0.8, 0.5), bone_liver = c(0.8, 0.5), lymph_liver = c(0, -0.85))
cols <- c(bone_lymph = "#2a78d6", bone_liver = "#eb6834", lymph_liver = "#1baf7a")
labs <- c(bone_lymph = "Bone vs Lymph node", bone_liver = "Bone vs Liver", lymph_liver = "Lymph node vs Liver")
pos <- data.table(region = c("unique_bone_lymph", "unique_bone_liver", "unique_lymph_liver", "site_bone",
                             "site_lymph_node", "site_liver", "shared_all_three"),
                  x = c(-1.6, 1.6, 0, 0, -0.95, 0.95, 0), y = c(1.0, 1.0, -1.85, 1.15, -0.45, -0.45, 0.05),
                  sub = c("", "", "", "bone-specific", "lymph-node-specific", "liver-specific", "all three"))
pos[, n := sapply(region, function(r) sum(mem$region == r))]
p <- ggplot()
for (k in names(cc)) p <- p + geom_polygon(data = circ(cc[[k]][1], cc[[k]][2]), aes(x, y), fill = cols[[k]],
                                           alpha = 0.18, colour = cols[[k]], linewidth = 0.8)
p <- p + geom_text(data = pos, aes(x, y, label = n), size = 6) +
    geom_text(data = pos[nzchar(sub)], aes(x, y - 0.28, label = sub), size = 3, colour = "#555555") +
    annotate("text", x = c(-1.9, 1.9, 0), y = c(2.25, 2.25, -2.6),
             label = sprintf("%s\n(%d genes)", labs, lengths(sig)[names(labs)]), size = 3.6,
             colour = cols[names(labs)]) +
    coord_equal(xlim = c(-3, 3), ylim = c(-3, 2.8)) + theme_void() +
    ggtitle(sprintf("Metastatic site differences: genes of significant stable regions (MethylDriver q < %g)", q_cut)) +
    theme(plot.title = element_text(size = 10, hjust = 0.5))
ggsave(file.path(out_dir, "venn_sites.pdf"), p, width = 7, height = 6.5)
ggsave(file.path(out_dir, "venn_sites.png"), p, width = 7, height = 6.5, dpi = 200)

# ---- GO BP ----
go_sets <- c("site_bone", "site_liver", "site_lymph_node", "shared_all_three", "all_site_genes")
tabs <- list()
for (s in go_sets) {
    g <- sets[[s]]
    if (is.null(g) || length(g) < min_genes) { cat(sprintf("[%s] %d genes - not tested\n", s, length(g))); next }
    ego <- tryCatch(enrichGO(g, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", universe = universe,
                             pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500),
                    error = function(e) NULL)
    df <- if (is.null(ego)) data.frame() else as.data.frame(ego)
    cat(sprintf("[%s] %d genes -> %d terms tested, %d with BH q < 0.05\n", s, length(g), nrow(df), sum(df$p.adjust < 0.05)))
    if (!nrow(df)) next
    fwrite(df, file.path(out_dir, sprintf("go_bp_%s.tsv", s)), sep = "\t")
    df$set <- s; tabs[[s]] <- df
    top <- head(df[order(df$pvalue), ], 15)
    top$Description <- factor(top$Description, levels = rev(top$Description))
    pp <- ggplot(top, aes(-log10(pvalue), Description, size = Count, colour = p.adjust < 0.05)) + geom_point() +
        scale_colour_manual(values = c(`TRUE` = "#2a78d6", `FALSE` = "#9a9a9a"), name = "BH q < 0.05") +
        labs(x = "-log10 p", y = NULL, title = sprintf("GO BP, %s (%d genes; universe %d)", s, length(g), length(universe))) +
        theme_minimal(base_size = 9)
    ggsave(file.path(out_dir, sprintf("go_bp_%s.pdf", s)), pp, width = 8, height = 5)
}
if (length(tabs)) fwrite(rbindlist(tabs, fill = TRUE), file.path(out_dir, "go_bp_all_sets.tsv"), sep = "\t")
