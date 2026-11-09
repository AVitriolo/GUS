# =============================================================================
# plot_demixt_vs_pames.R
# DeMixT RNA tumour proportion (pi) vs PAMES methylation purity, per CPGEA
# Tumour. The two measure different things (share of RNA vs share of DNA),
# so they should correlate but need not match; large disagreements flag
# Tumours whose purified RNA and purified methylation describe different
# tumour content.
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel) })

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

pi_paths <- strsplit(args$pi, ",")[[1]]      # demixt_pi_<split>.tsv, one or more
pames    <- args$pames
out_png  <- args$out_png
n_label  <- as.integer(if (is.null(args$n_label)) 8 else args$n_label)

pi <- rbindlist(lapply(pi_paths, function(f)
    fread(f)[, split := sub("^testjoint$", "test", sub("_ns.*", "", sub("^demixt_pi_(.*)\\.tsv$", "\\1", basename(f))))]))
pm <- fread(pames)
setnames(pm, c("sample", "pames", "condition"))
d  <- merge(pi, pm[, .(sample, pames)], by = "sample")
stopifnot(nrow(d) == nrow(pi))

d[, diff := pi_tumour - pames]
lab <- d[order(-abs(diff))][seq_len(min(n_label, .N))]

r_s <- cor(d$pi_tumour, d$pames, method = "spearman")
r_p <- cor(d$pi_tumour, d$pames)
cat(sprintf("n = %d | Spearman %.2f | Pearson %.2f | median DeMixT %.2f vs PAMES %.2f | median |diff| %.2f\n",
            nrow(d), r_s, r_p, median(d$pi_tumour), median(d$pames), median(abs(d$diff))))
print(lab[, .(sample, split, pi_tumour = round(pi_tumour, 2), pames = round(pames, 2), diff = round(diff, 2))])

INK <- "#3d3d3a"; MUTED <- "#8a8984"; GRID <- "#e6e5df"
p <- ggplot(d, aes(pames, pi_tumour)) +
    geom_abline(slope = 1, intercept = 0, colour = MUTED, linewidth = 0.4, linetype = "dashed") +
    geom_point(aes(shape = split), colour = "#2a78d6", size = 2.4, alpha = 0.85) +
    geom_text_repel(data = lab, aes(label = sample), colour = INK, size = 3, min.segment.length = 0, segment.colour = MUTED, box.padding = 0.5, seed = 1) +
    scale_shape_manual(values = c(train = 16, test = 17), name = NULL) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "PAMES purity (methylation, DNA)", y = "DeMixT tumour proportion (RNA)",
         title = "Tumour purity: DeMixT (RNA) vs PAMES (methylation)",
         subtitle = sprintf("%d CPGEA Tumours. Spearman %.2f, Pearson %.2f\nDashed = equal. Labelled: the %d largest disagreements",
                            nrow(d), r_s, r_p, nrow(lab))) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
          axis.text = element_text(colour = MUTED), axis.title = element_text(colour = INK),
          plot.title = element_text(colour = INK, face = "bold"),
          plot.subtitle = element_text(colour = MUTED),
          legend.position = if (length(unique(d$split)) > 1) "bottom" else "none",
          plot.background = element_rect(fill = "white", colour = NA))

ggsave(out_png, p, width = 7, height = 7.4, dpi = 150)
ggsave(sub("\\.png$", ".pdf", out_png), p, width = 7, height = 7.4, device = cairo_pdf)
cat("->", out_png, "\n")
