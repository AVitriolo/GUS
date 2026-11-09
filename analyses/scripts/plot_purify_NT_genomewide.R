# =============================================================================
# plot_purify_NT_genomewide.R
# Genome-wide view of the CAMDAC-equation purification of CPGEA Tumours
# (purify_tumour_NT.R), on a random sample of CpG rows:
#   page 1: bulk vs purified beta (hexbin), inside called peaks vs outside
#   page 2: distribution of purified - bulk, inside vs outside
#   page 3: per-Tumour purity vs fraction clipped / mean |change| (all CpGs)
# =============================================================================
suppressPackageStartupMessages({ library(GenomicRanges); library(rhdf5); library(ggplot2) })

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

orig_h5   <- args$orig_h5
pur_h5    <- args$pur_h5
rowranges <- args$rowranges
coldata   <- args$coldata
peaks     <- args$peaks
summary_p <- args$summary
n_rows    <- as.integer(args$n_rows)
out_pdf   <- args$out_pdf

set.seed(1)
rr <- readRDS(rowranges)
cd <- as.data.frame(readRDS(coldata))
pk <- read.delim(peaks)
sm <- read.delim(summary_p)

t_col <- which(cd$condition == "Tumour")
# random blocks of consecutive CpGs: few hyperslab reads instead of n_rows
BLOCK  <- 1000L
starts <- sort(sample(seq(1L, length(rr) - BLOCK + 1L, by = BLOCK), n_rows %/% BLOCK))
rows   <- as.vector(outer(seq_len(BLOCK) - 1L, starts, "+"))

# contiguous runs of rows so h5read gets few hyperslabs
read_rows <- function(h5, rows, cols) {
    runs <- split(rows, cumsum(c(1, diff(rows) != 1)))
    do.call(rbind, lapply(runs, function(r) h5read(h5, "beta", index = list(r, cols))))
}
bulk <- read_rows(orig_h5, rows, t_col)
pur  <- read_rows(pur_h5,  rows, t_col)

pk   <- pk[pk$sample %in% rownames(cd)[t_col], ]
hits <- findOverlaps(rr[rows], GRanges(pk$seqnames, IRanges(pk$start, pk$end)))
in_peak <- matrix(FALSE, nrow(bulk), ncol(bulk))
in_peak[cbind(queryHits(hits), match(pk$sample[subjectHits(hits)], rownames(cd)[t_col]))] <- TRUE

df <- data.frame(bulk = as.vector(bulk), purified = as.vector(pur),
                 panel = ifelse(as.vector(in_peak), "Inside called peaks (n_t from CN)",
                                "Outside peaks (n_t = 2)"))
df$panel <- factor(df$panel, levels = c("Outside peaks (n_t = 2)", "Inside called peaks (n_t from CN)"))
df$delta <- df$purified - df$bulk

stats <- do.call(rbind, lapply(split(df, df$panel), function(d)
    data.frame(panel = d$panel[1], n = nrow(d),
               med = median(abs(d$delta)), p90 = quantile(abs(d$delta), 0.9),
               clip = mean(d$purified %in% c(0, 1) & !(d$bulk %in% c(0, 1))))))
stats$label <- sprintf("n = %s cells\nmedian |Δβ| = %.3f\n90th pct |Δβ| = %.3f\nclipped to 0/1 = %.1f%%",
                       format(stats$n, big.mark = ","), stats$med, stats$p90, 100 * stats$clip)
print(stats[, c("panel", "n", "med", "p90", "clip")], row.names = FALSE)

INK   <- "#3d3d3a"; MUTED <- "#8a8984"; GRID <- "#e6e5df"
theme_set(theme_minimal(base_size = 11) +
          theme(panel.grid.minor = element_blank(),
                panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
                strip.text = element_text(face = "bold", colour = INK, hjust = 0),
                axis.text = element_text(colour = MUTED),
                axis.title = element_text(colour = INK),
                plot.title = element_text(colour = INK, face = "bold"),
                plot.subtitle = element_text(colour = MUTED),
                plot.background = element_rect(fill = "white", colour = NA)))
sub <- sprintf("%s CpGs (%d random blocks of %d) x %d Tumours, PAMES purity, matched Normal as b_n",
               format(length(rows), big.mark = ","), length(starts), BLOCK, length(t_col))

p1 <- ggplot(df, aes(bulk, purified)) +
    geom_abline(slope = 1, intercept = 0, colour = MUTED, linewidth = 0.4, linetype = "dashed") +
    geom_hex(bins = 70) +
    scale_fill_gradient(low = "#cde2fb", high = "#0d366b", trans = "log10", name = "cells") +
    facet_wrap(~ panel, ncol = 2) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "Bulk beta (original)", y = "Purified tumour beta",
         title = "Genome-wide: bulk vs purified Tumour beta", subtitle = paste0(sub, ". Dashed = no change"))

p2 <- ggplot(df, aes(delta)) +
    geom_vline(xintercept = 0, colour = MUTED, linewidth = 0.4) +
    geom_histogram(binwidth = 0.01, fill = "#2a78d6", colour = "white", linewidth = 0.1) +
    geom_text(data = stats, aes(x = -Inf, y = Inf, label = label), hjust = -0.05, vjust = 1.2,
              size = 3, colour = INK, lineheight = 0.95) +
    facet_wrap(~ panel, ncol = 2, scales = "free_y") +
    coord_cartesian(xlim = c(-0.6, 0.6)) +
    labs(x = "Δβ = purified - bulk", y = "Tumour CpG × sample cells",
         title = "Genome-wide: how much the Tumour betas changed", subtitle = paste0(sub, ". x-axis limited to ±0.6"))

sm_long <- rbind(data.frame(rho = sm$rho, value = 100 * sm$frac_clipped, metric = "Clipped to 0/1 (%)"),
                 data.frame(rho = sm$rho, value = sm$mean_abs_delta,    metric = "Mean |purified - bulk|"))
p3 <- ggplot(sm_long, aes(rho, value)) +
    geom_point(colour = "#2a78d6", size = 2.2, alpha = 0.85) +
    facet_wrap(~ metric, ncol = 2, scales = "free_y") +
    labs(x = "PAMES purity", y = NULL,
         title = "Per Tumour (all 19M CpGs): lower purity -> larger change and more clipping",
         subtitle = sprintf("%d Tumours. Spearman: clipped %.2f, mean change %.2f", nrow(sm),
                            cor(sm$rho, sm$frac_clipped, method = "s"),
                            cor(sm$rho, sm$mean_abs_delta, method = "s")))

cairo_pdf(out_pdf, onefile = TRUE, width = 10, height = 6)
print(p1); print(p2); print(p3)
invisible(dev.off())
ggsave(sub("\\.pdf$", "_scatter.png",    out_pdf), p1, width = 10, height = 6, dpi = 150)
ggsave(sub("\\.pdf$", "_delta.png",      out_pdf), p2, width = 10, height = 5, dpi = 150)
ggsave(sub("\\.pdf$", "_per_tumour.png", out_pdf), p3, width = 10, height = 4.5, dpi = 150)
cat("->", out_pdf, "\n")
