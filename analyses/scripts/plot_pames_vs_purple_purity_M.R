# =============================================================================
# plot_pames_vs_purple_purity_M.R
# PAMES (WGBS) vs PURPLE (WGS, Lundberg et al. WCDT) purity for the metastases.
# PURPLE is read from the per-site table already extracted from the poppy RData.
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2) })

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
pames_path  <- args$pames_table
purple_path <- args$purple_table
out_pdf     <- args$out_pdf

p <- read.delim(pames_path)
w <- read.delim(purple_path)
names(w)[names(w) == "purity"] <- "purple"
d <- merge(p[, c("sample", "pames", "site")], w[, c("sample", "purple")], by = "sample")

r   <- cor(d$pames, d$purple)
rho <- cor(d$pames, d$purple, method = "spearman")
cat(sprintf("n = %d | Pearson r = %.3f | Spearman rho = %.3f\n", nrow(d), r, rho))

INK <- "#3d3d3a"; MUTED <- "#8a8984"; GRID <- "#e6e5df"
th <- theme_minimal(base_size = 8, base_family = "DejaVu Sans") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
          axis.text = element_text(colour = MUTED), axis.title = element_text(colour = INK),
          plot.title = element_text(colour = INK), plot.subtitle = element_text(colour = MUTED),
          plot.background = element_rect(fill = "white", colour = NA))

g <- ggplot(d, aes(purple, pames)) +
    geom_smooth(method = "lm", se = TRUE, colour = "#2a78d6", linewidth = 0.6) +
    geom_point(colour = INK, size = 1.3, alpha = 0.75) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "WGS purity (PURPLE, Lundberg et al. WCDT)", y = "WGBS purity (PAMES)",
         title = "Metastatic (DTB) sample purity: PAMES vs WGS-PURPLE",
         subtitle = sprintf("n = %d, Pearson r = %.2f, Spearman rho = %.2f", nrow(d), r, rho)) + th

cairo_pdf(out_pdf, width = 4.6, height = 4.2); print(g); invisible(dev.off())
ggsave(sub("\\.pdf$", ".png", out_pdf), g, width = 4.6, height = 4.2, dpi = 300)
cat("->", out_pdf, "\n")
