#!/usr/bin/env Rscript
# =============================================================================
# plot_cpg_annotation.R
# Plots of annotate_selected_cpgs.py: % of selected CpGs vs % of candidate CpGs
#   genomic_annotation.png : promoter / proximal / distal enhancer (ENCODE cCREs),
#                            exon / intron / intergenic (one category per CpG)
#   repeat_annotation.png  : RepeatMasker families grouped by class (families
#                            < --min_pct of the candidates folded into "other")
# Ratio selected / candidates written next to each pair of bars.
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
args    <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
in_dir  <- args$in_dir
min_pct <- as.numeric(if (is.null(args$min_pct)) 0.3 else args$min_pct)

cols <- c("Selected CpGs" = "#4a4a48", "Candidate CpGs" = "#c9c8c3")    # dark / light grey
ink  <- c(primary = "#0b0b0b", secondary = "#52514e", grid = "#e6e5e0")
set_lab <- c(selected = "Selected CpGs", candidates = "Candidate CpGs")
theme_ann <- theme_minimal(base_size = 10) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          panel.grid.major.x = element_line(colour = ink[["grid"]], linewidth = 0.3),
          axis.text.y = element_text(colour = ink[["primary"]]), axis.text.x = element_text(colour = ink[["secondary"]]),
          axis.title.x = element_text(colour = ink[["secondary"]], size = 9),
          legend.position = "top", legend.justification = "left", legend.title = element_blank(),
          plot.background = element_rect(fill = "white", colour = NA))

pair_plot <- function(d, ylab_order, facet = NULL) {
    d <- copy(d)[set %in% names(set_lab)]
    d[, set := factor(set_lab[set], levels = rev(set_lab))]
    d[, label := factor(label, levels = rev(ylab_order))]
    r <- dcast(d, as.formula(paste(paste(c("label", facet), collapse = "+"), "~ set")), value.var = "pct")
    r[, ratio := `Selected CpGs` / `Candidate CpGs`]
    r[, xpos := pmax(`Selected CpGs`, `Candidate CpGs`) + 0.012 * max(d$pct)]
    p <- ggplot(d, aes(x = pct, y = label, fill = set)) +
        geom_col(position = position_dodge(width = 0.8), width = 0.75) +
        geom_text(data = r, aes(x = xpos, y = label, label = sprintf("%.2f×", ratio)), inherit.aes = FALSE,
                  hjust = 0, size = 2.9, colour = ink[["secondary"]]) +
        scale_fill_manual(values = cols, breaks = unname(set_lab)) +
        scale_x_continuous(expand = expansion(mult = c(0, 0.14))) +
        labs(x = "% of CpGs", y = NULL) + theme_ann
    if (!is.null(facet)) p <- p + facet_grid(as.formula(paste(facet, "~ .")), scales = "free_y", space = "free_y", switch = "y",
                                             labeller = as_labeller(function(x) ifelse(x %in% c("No repeat", "Other"), "", x))) +
        theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", colour = ink[["primary"]]),
              panel.spacing.y = unit(6, "pt"))
    p
}

# ---- genomic annotation ----
g <- fread(file.path(in_dir, "genomic_annotation_counts.tsv"))
g[, label := category]
gl <- c("Promoter (PLS)", "Proximal enhancer (pELS)", "Distal enhancer (dELS)", "Exon", "Intron", "Intergenic")
p1 <- pair_plot(g, gl)
ggsave(file.path(in_dir, "genomic_annotation.png"), p1, width = 7, height = 3.6, dpi = 250, bg = "white")
ggsave(file.path(in_dir, "genomic_annotation.pdf"), p1, width = 7, height = 3.6, device = cairo_pdf)

# ---- repeats: families grouped by class ----
f <- fread(file.path(in_dir, "repeat_family_counts.tsv"))
f <- f[set %in% names(set_lab)]
f[, repClass := sub("\\?$", "", repClass)]; f[, repFamily := sub("\\?$", "", repFamily)]
f[, cls := fcase(repClass %in% c("LINE", "SINE", "LTR"), repClass, repClass == "DNA", "DNA transposon",
                 repClass %in% c("Simple_repeat", "Low_complexity", "Satellite"), "Simple / low\ncomplexity / satellite",
                 repClass == "No repeat", "No repeat", default = "Other")]
big <- f[set == "candidates", .(pct = sum(pct)), by = .(cls, repFamily)][pct >= min_pct & !cls %in% c("Other", "No repeat")]
f[, fam := fifelse(paste(cls, repFamily) %in% big[, paste(cls, repFamily)], repFamily,
           fifelse(cls %in% c("Other", "No repeat"), cls, paste("other", cls)))]
f[cls == "Simple / low\ncomplexity / satellite", fam := sub("_", " ", repClass)]
r <- f[, .(pct = sum(pct)), by = .(set, cls, label = fam)]
cls_order <- c("No repeat", "LINE", "SINE", "LTR", "DNA transposon", "Simple / low\ncomplexity / satellite", "Other")
r[, cls := factor(cls, levels = cls_order)]
fam_order <- r[set == "candidates"][order(cls, grepl("^other", label), -pct)][, unique(label)]
tot <- f[repClass != "No repeat", .(pct = sum(pct)), by = set]
p2 <- pair_plot(r, fam_order, facet = "cls") +
    labs(caption = sprintf("In a repeat: %.1f%% of selected vs %.1f%% of candidate CpGs",
                           tot[set == "selected", pct], tot[set == "candidates", pct])) +
    theme(plot.caption = element_text(colour = ink[["secondary"]], hjust = 0, size = 8.5))

# ---- pie on the right: genomic regions of the selected CpGs that lie in a repeat ----
a <- fread(file.path(in_dir, "selected_cpgs_annotation.tsv.gz"), select = c("category", "repClass"))
pie_d <- a[repClass != "No repeat", .N, by = category][, category := factor(category, levels = gl)][order(category)]
pie_d[, pct := 100 * N / sum(N)]
pie_d[, lab := factor(sprintf("%s  %.1f%%", category, pct), levels = sprintf("%s  %.1f%%", category, pct))]
ramp <- c("#0d366b", "#1c5cab", "#3987e5", "#6da7ec", "#9ec5f4", "#cde2fb")      # promoter (dark) -> intergenic (light)
pie <- ggplot(pie_d, aes(x = 1, y = N, fill = lab)) +
    geom_col(width = 1, colour = "white", linewidth = 0.4) +
    coord_polar(theta = "y", direction = -1) +
    scale_fill_manual(values = setNames(ramp, levels(pie_d$lab)), name = NULL) +
    labs(title = sprintf("Genomic region of the selected\nCpGs in repeats (n = %s)", format(sum(pie_d$N), big.mark = ","))) +
    theme_void(base_size = 9) +
    theme(legend.position = "bottom", legend.direction = "vertical", legend.key.size = unit(10, "pt"),
          legend.text = element_text(size = 8, colour = ink[["primary"]]),
          plot.title = element_text(size = 9, hjust = 0.5, colour = ink[["primary"]]),
          plot.background = element_rect(fill = "white", colour = NA))

draw_combo <- function() {
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2, widths = grid::unit(c(7.5, 3), "in"))))
    print(p2, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
    grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
    print(pie, vp = grid::viewport(x = 0.5, y = 0.62, width = 0.95, height = 0.62))
    grid::popViewport(2)
}
png(file.path(in_dir, "repeat_annotation.png"), width = 10.5, height = 6.2, units = "in", res = 250, bg = "white", type = "cairo")
draw_combo(); invisible(dev.off())
cairo_pdf(file.path(in_dir, "repeat_annotation.pdf"), width = 10.5, height = 6.2)
draw_combo(); invisible(dev.off())

# ---- class-level version (LINE, SINE, LTR, DNA, ...), same pie ----
rc <- f[, .(pct = sum(pct)), by = .(set, label = gsub("\n", " ", cls))]
p3 <- pair_plot(rc, gsub("\n", " ", cls_order)) +
    labs(caption = sprintf("In a repeat: %.1f%% of selected vs %.1f%% of candidate CpGs",
                           tot[set == "selected", pct], tot[set == "candidates", pct])) +
    theme(plot.caption = element_text(colour = ink[["secondary"]], hjust = 0, size = 8.5))
draw_combo_cls <- function() {
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2, widths = grid::unit(c(7, 3), "in"))))
    print(p3, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
    grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
    print(pie, vp = grid::viewport(x = 0.5, y = 0.5, width = 0.95, height = 0.95))
    grid::popViewport(2)
}
png(file.path(in_dir, "repeat_annotation_classes.png"), width = 10, height = 4.2, units = "in", res = 250, bg = "white", type = "cairo")
draw_combo_cls(); invisible(dev.off())
cairo_pdf(file.path(in_dir, "repeat_annotation_classes.pdf"), width = 10, height = 4.2)
draw_combo_cls(); invisible(dev.off())
cat("plots written to", in_dir, "\n")
