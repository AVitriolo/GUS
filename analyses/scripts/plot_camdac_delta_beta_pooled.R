# =============================================================================
# plot_camdac_delta_beta_pooled.R
# Same delta as plot_camdac_delta_beta.R (m_t purified - m bulk), but pooled
# over all CpGs instead of split into diploid vs CNA peaks. Reads the histogram
# counts that plot_camdac_delta_beta.R already wrote, so the pure.csv.gz files
# are never re-read.
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
hist_tsv <- args$histogram_tsv
out_pdf  <- args$out_pdf
label    <- if (is.null(args$label)) "Tumours" else args$label

h <- fread(hist_tsv)
h[, total := diploid + cna]
n_tot <- sum(as.numeric(h$total))

# median |delta|, read off the pooled counts (bins ordered by distance from 0)
w <- h[order(abs(bin))]
med_abs <- abs(w$bin[which(cumsum(as.numeric(w$total)) >= n_tot / 2)[1]])
within <- function(x) 100 * sum(as.numeric(h[abs(bin) <= x, total])) / n_tot

h[, share := as.numeric(total) / n_tot]

cat(sprintf("%s CpG x Tumour values | median |delta| ~%.3f | within +-0.05: %.1f%% | within +-0.1: %.1f%%\n",
            format(n_tot, big.mark = ","), med_abs, within(0.05), within(0.1)))

INK <- "#3d3d3a"; MUTED <- "#8a8984"; GRID <- "#e6e5df"
th <- theme_minimal(base_size = 8, base_family = "DejaVu Sans") +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
          axis.text = element_text(colour = MUTED), axis.title = element_text(colour = INK),
          plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = MUTED),
          plot.background = element_rect(fill = "white", colour = NA))

p <- ggplot(h, aes(bin, share)) +
    geom_vline(xintercept = 0, colour = MUTED, linewidth = 0.4) +
    geom_col(fill = "#2a78d6", width = 0.01) +
    scale_y_continuous(labels = function(x) paste0(round(100 * x, 1), "%")) +
    coord_cartesian(xlim = c(-0.6, 0.6)) +
    labs(x = "Δβ = CAMDAC purified tumour beta - bulk beta",
         y = "% of CpGs (all tumoral samples)") + th

cairo_pdf(out_pdf, width = 4.2, height = 3.2); print(p); invisible(dev.off())
ggsave(sub("\\.pdf$", ".png", out_pdf), p, width = 4.2, height = 3.2, dpi = 300)
cat("->", out_pdf, "\n")
