# =============================================================================
# plot_cna_correction_test.R
# How much the Tumour betas change on the small test h5 (rows around three
# peaks: chr8 amp peak 9, chr16 del peak 16, chr3 del peak 4), PAMES purity:
#   - CNA-only correction  (camdac_correct_NT.R, re-mixed at n_t = 2)
#   - purified tumour beta (purify_tumour_NT.R)
# Page 1: bulk vs corrected beta. Page 2: distribution of corrected - bulk.
# =============================================================================
suppressPackageStartupMessages({ library(GenomicRanges); library(rhdf5); library(ggplot2) })

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

orig_h5   <- args$orig_h5
cna_h5    <- args$cna_h5
pur_h5    <- args$pur_h5
rowranges <- args$rowranges
coldata   <- args$coldata
peaks     <- args$peaks
out_pdf   <- args$out_pdf

rr <- readRDS(rowranges)
cd <- as.data.frame(readRDS(coldata))
pk <- read.delim(peaks)

bulk <- h5read(orig_h5, "beta")
cna  <- h5read(cna_h5,  "beta")
pur  <- h5read(pur_h5,  "beta")
t_col <- which(cd$condition == "Tumour")

# (row, Tumour column) cells inside a peak called in that sample
pk   <- pk[pk$sample %in% rownames(cd)[t_col], ]
hits <- findOverlaps(rr, GRanges(pk$seqnames, IRanges(pk$start, pk$end)))
called <- data.frame(row = queryHits(hits),
                     col = match(pk$sample[subjectHits(hits)], rownames(cd)),
                     type = pk$type[subjectHits(hits)])
in_peak <- matrix(FALSE, nrow(bulk), ncol(bulk))
in_peak[cbind(called$row, called$col)] <- TRUE

cells <- function(idx, corrected, label)
    data.frame(bulk = bulk[idx], corrected = corrected[idx], panel = label)

all_t <- as.matrix(expand.grid(row = seq_len(nrow(bulk)), col = t_col))
df <- rbind(
    cells(as.matrix(called[called$type == "amp", 1:2]), cna, "CNA-only: amplification peaks"),
    cells(as.matrix(called[called$type == "del", 1:2]), cna, "CNA-only: deletion peaks"),
    cells(all_t[ in_peak[all_t], , drop = FALSE],       pur, "Purified: inside called peaks"),
    cells(all_t[!in_peak[all_t], , drop = FALSE],       pur, "Purified: outside peaks")
)
df$panel <- factor(df$panel, levels = unique(df$panel))
df$delta <- df$corrected - df$bulk

stats <- do.call(rbind, lapply(split(df, df$panel), function(d)
    data.frame(panel = d$panel[1], n = nrow(d),
               med = median(abs(d$delta)), p90 = quantile(abs(d$delta), 0.9),
               clip = mean(d$corrected %in% c(0, 1) & !(d$bulk %in% c(0, 1))))))
stats$label <- sprintf("n = %s cells\nmedian |Δβ| = %.3f\n90th pct |Δβ| = %.3f",
                       format(stats$n, big.mark = ","), stats$med, stats$p90)
# CNA-only re-mixes after clipping, so a clipped cell is not visible in its
# output; the clipped share only means something for the purified beta
is_pur <- grepl("^Purified", stats$panel)
stats$label[is_pur] <- sprintf("%s\nclipped to 0/1 = %.1f%%", stats$label[is_pur], 100 * stats$clip[is_pur])
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

p1 <- ggplot(df, aes(bulk, corrected)) +
    geom_abline(slope = 1, intercept = 0, colour = MUTED, linewidth = 0.4, linetype = "dashed") +
    geom_hex(bins = 60) +
    scale_fill_gradient(low = "#cde2fb", high = "#0d366b", trans = "log10", name = "cells") +
    facet_wrap(~ panel, ncol = 2) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "Bulk beta (original)", y = "Corrected beta",
         title = "Tumour beta before vs after correction (test h5, PAMES purity)",
         subtitle = "Dashed line = no change. Rows around chr8 amp peak 9, chr16 del peak 16, chr3 del peak 4")

p2 <- ggplot(df, aes(delta)) +
    geom_vline(xintercept = 0, colour = MUTED, linewidth = 0.4) +
    geom_histogram(binwidth = 0.01, fill = "#2a78d6", colour = "white", linewidth = 0.1) +
    geom_text(data = stats, aes(x = -Inf, y = Inf, label = label), hjust = -0.05, vjust = 1.2,
              size = 3, colour = INK, lineheight = 0.95) +
    facet_wrap(~ panel, ncol = 2, scales = "free_y") +
    coord_cartesian(xlim = c(-0.5, 0.5)) +
    labs(x = "Δβ = corrected - bulk", y = "Tumour CpG × sample cells",
         title = "How much the Tumour betas changed (test h5, PAMES purity)",
         subtitle = "x-axis limited to ±0.5")

cairo_pdf(out_pdf, onefile = TRUE, width = 9, height = 8)
print(p1); print(p2)
invisible(dev.off())
ggsave(sub("\\.pdf$", "_scatter.png", out_pdf), p1, width = 9, height = 8, dpi = 150)
ggsave(sub("\\.pdf$", "_delta.png",   out_pdf), p2, width = 9, height = 7, dpi = 150)
cat("->", out_pdf, "\n")
