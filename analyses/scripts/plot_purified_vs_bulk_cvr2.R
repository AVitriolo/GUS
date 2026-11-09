# =============================================================================
# plot_purified_vs_bulk_cvr2.R
# Paired per-transcript cv_r2: models on purity-corrected inputs (CAMDAC betas
# + bMIND expression) vs the uncorrected kcapped baseline.
#   left : joint density of the two cv_r2 values, line of equality
#   right: distribution of the per-transcript difference
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
m <- read.delim(args$paired_tsv)
out_pdf <- args$out_pdf

d <- m$purified - m$bulk
stat <- sprintf("n = %s transcripts   median Δ = %+.4f   improved %.0f%% / worsened %.0f%%",
                format(nrow(m), big.mark = ","), median(d), 100*mean(d>0), 100*mean(d<0))
cat(stat, "\n")

INK <- "#3d3d3a"; MUTED <- "#8a8984"; GRID <- "#e6e5df"
th <- theme_minimal(base_size = 8, base_family = "DejaVu Sans") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
          axis.text = element_text(colour = MUTED), axis.title = element_text(colour = INK),
          plot.title = element_text(colour = INK), plot.subtitle = element_text(colour = MUTED),
          legend.title = element_text(colour = INK), legend.text = element_text(colour = MUTED),
          plot.background = element_rect(fill = "white", colour = NA))

lim <- quantile(c(m$purified, m$bulk), c(0.001, 0.999))
p1 <- ggplot(m, aes(bulk, purified)) +
    geom_bin2d(bins = 90) +
    scale_fill_gradient(low = "#dce8f7", high = "#1a5fb4", trans = "log10", name = "transcripts") +
    geom_abline(slope = 1, intercept = 0, colour = INK, linewidth = 0.4) +
    coord_cartesian(xlim = lim, ylim = lim) +
    labs(x = "cv R², uncorrected (bulk)", y = "cv R², purity-corrected",
         subtitle = sprintf("Spearman %.3f; line = equal performance",
                            cor(m$purified, m$bulk, method = "spearman"))) + th

p2 <- ggplot(data.frame(d = d), aes(d)) +
    geom_histogram(bins = 90, fill = "#2a78d6", colour = NA) +
    geom_vline(xintercept = 0, colour = INK, linewidth = 0.4) +
    coord_cartesian(xlim = quantile(d, c(0.005, 0.995))) +
    labs(x = "Δ cv R²  (corrected − uncorrected)", y = "transcripts",
         subtitle = sprintf("median %+.4f; middle 50%% within %+.3f to %+.3f",
                            median(d), quantile(d,.25), quantile(d,.75))) + th

g <- (p1 | p2) + patchwork::plot_annotation(
        title = "Purity correction does not change prediction accuracy",
        subtitle = stat,
        theme = th)
cairo_pdf(out_pdf, width = 8.6, height = 3.6); print(g); invisible(dev.off())
ggsave(sub("\\.pdf$", ".png", out_pdf), g, width = 8.6, height = 3.6, dpi = 300)
cat("->", out_pdf, "\n")
