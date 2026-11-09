# =============================================================================
# plot_camdac_delta_beta_per_tumour.R
# Per CPGEA Tumour vs PAMES purity: median |delta beta| and mean delta beta
# (delta = CAMDAC purified m_t - bulk m), from the per-Tumour table written by
# plot_camdac_delta_beta.R (no need to re-read the pure.csv.gz files).
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
stats   <- fread(args$stats)
out_png <- args$out_png

sl <- rbind(data.table(purity = stats$purity, value = stats$median_absdelta, metric = "Median |Δβ| per Tumour"),
            data.table(purity = stats$purity, value = stats$mean_delta,      metric = "Mean Δβ per Tumour"))
sl[, metric := factor(metric, levels = c("Median |Δβ| per Tumour", "Mean Δβ per Tumour"))]
metrics <- if (is.null(args$metrics)) levels(sl$metric) else strsplit(args$metrics, ",")[[1]]
sl <- sl[metric %in% metrics]
sl[, metric := droplevels(metric)]
rho <- sl[, .(r = cor(purity, value, method = "spearman")), by = metric]
rho[, label := sprintf("Spearman with purity %.2f", r)]

INK <- "#3d3d3a"; MUTED <- "#8a8984"; GRID <- "#e6e5df"
p <- ggplot(sl, aes(purity, value)) +
    { if ("Mean Δβ per Tumour" %in% sl$metric)
          geom_hline(data = data.table(metric = factor("Mean Δβ per Tumour", levels = levels(sl$metric)), y = 0),
                     aes(yintercept = y), colour = MUTED, linewidth = 0.4) } +
    geom_point(colour = "#2a78d6", size = 2.2, alpha = 0.85) +
    geom_text(data = rho, aes(x = -Inf, y = Inf, label = label), hjust = -0.05, vjust = 1.5,
              size = 3.2, colour = INK) +
    facet_wrap(~ metric, ncol = nlevels(sl$metric), scales = "free_y") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +   # headroom for the label
    labs(x = "PAMES purity", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
          strip.text = element_text(face = "bold", colour = INK, hjust = 0),
          axis.text = element_text(colour = MUTED), axis.title = element_text(colour = INK),
          plot.background = element_rect(fill = "white", colour = NA))

w <- if (nlevels(sl$metric) == 1) 6 else 11
ggsave(out_png, p, width = w, height = 4.4, dpi = 150)
ggsave(sub("\\.png$", ".pdf", out_png), p, width = w, height = 4.4, device = cairo_pdf)
cat("->", out_png, "\n")
