# =============================================================================
# plot_demixt_delta_expression.R
# How much DeMixT changed the CPGEA Tumour expression:
#   log2FC = log2((purified + 1) / (bulk + 1)), per transcript x Tumour, with
#   bulk = the upper-quartile-scaled (+1) counts DeMixT was given.
#   page 1: pooled log2FC distribution
#   page 2: per Tumour median |log2FC| vs DeMixT tumour proportion
#   page 3: marker genes (transcripts summed per gene): non-tumour genes
#           (stroma / smooth muscle / basal / immune) should go down, prostate
#           tumour genes up
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
dir     <- args$dir
tags    <- strsplit(args$tags, ",")[[1]]       # e.g. train_ns10_ng1500,testjoint_ns14_ng1500
conv    <- args$conversion
out_pdf <- args$out_pdf

MARKERS <- data.table(
    gene  = c("MYH11", "ACTA2", "DES", "CNN1", "COL1A1", "COL3A1", "KRT5", "KRT15", "TP63", "PTPRC",
              "AMACR", "PCA3", "HPN", "FOLH1", "DLX1", "TDRD1", "ERG"),
    group = c(rep("Non-tumour: expected down", 10), rep("Prostate tumour: expected up", 7)))

load_tag <- function(tag) {
    pur  <- fread(file.path(dir, sprintf("demixt_purified_tumour_%s.tsv.gz", tag)))
    bulk <- fread(file.path(dir, sprintf("bulk_scaled_plus1_tumour_%s.tsv.gz", tag)))
    pi   <- fread(file.path(dir, sprintf("demixt_pi_%s.tsv", tag)))
    bulk <- bulk[match(pur$transcript_id, bulk$transcript_id)]
    stopifnot(identical(names(pur), names(bulk)))
    list(tx = pur$transcript_id, P = as.matrix(pur[, -1]), B = as.matrix(bulk[, -1]),
         pi = pi, split = sub("^testjoint$", "test", sub("_(ns.*|givenpi)$", "", tag)))
}
dat <- lapply(tags, load_tag)

# ---- pooled distribution + per Tumour ----
BREAKS <- seq(-8, 8, by = 0.1)
hist <- rbindlist(lapply(dat, function(d) {
    fc <- log2((d$P + 1) / (d$B + 1))
    data.table(split = d$split, bin = head(BREAKS, -1) + 0.05,
               n = as.numeric(table(cut(pmin(pmax(fc, -7.99), 7.99), BREAKS))))
}))
per_t <- rbindlist(lapply(dat, function(d) {
    fc <- log2((d$P + 1) / (d$B + 1))
    data.table(split = d$split, sample = colnames(d$P),
               pi = d$pi$pi_tumour[match(colnames(d$P), d$pi$sample)],
               median_abs_log2fc = apply(abs(fc), 2, median),
               frac_down = colMeans(fc < 0))
}))
fwrite(per_t, sub("\\.pdf$", "_per_tumour.tsv", out_pdf), sep = "\t")

# ---- marker genes: sum transcripts per gene, then log2FC ----
cv <- fread(conv)[, .(transcript_id, gene_name)]
mk <- rbindlist(lapply(dat, function(d) {
    g <- cv$gene_name[match(d$tx, cv$transcript_id)]
    keep <- g %in% MARKERS$gene
    if (!any(keep)) return(NULL)
    Pg <- rowsum(d$P[keep, , drop = FALSE], g[keep]); Bg <- rowsum(d$B[keep, , drop = FALSE], g[keep])
    fc <- log2((Pg + 1) / (Bg + 1))
    melt(data.table(gene = rownames(fc), fc), id.vars = "gene", variable.name = "sample", value.name = "log2fc")[, split := d$split]
}))
mk <- merge(mk, MARKERS, by = "gene")
mk[, gene := factor(gene, levels = MARKERS$gene)]
missing <- setdiff(MARKERS$gene, unique(as.character(mk$gene)))

pool <- hist[, .(n = sum(n)), by = bin]
w <- pool[order(abs(bin))]
cat(sprintf("%s transcript x Tumour values | %.1f%% decreased | median |log2FC| ~%.2f\n",
            format(sum(pool$n), big.mark = ","), 100 * sum(pool[bin < 0, n]) / sum(pool$n),
            abs(w$bin[which(cumsum(w$n) >= sum(w$n) / 2)[1]])))
mk[, expected := if (group[1] == "Non-tumour: expected down") log2fc < 0 else log2fc > 0, by = group]
print(mk[, .(median_log2fc = round(median(log2fc), 2), frac_in_expected_direction = round(mean(expected), 2)),
         by = .(group, gene)][order(group, gene)])
if (length(missing)) cat("markers not among the deconvolved transcripts:", paste(missing, collapse = ", "), "\n")

INK <- "#3d3d3a"; MUTED <- "#8a8984"; GRID <- "#e6e5df"
th <- theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
          strip.text = element_text(face = "bold", colour = INK, hjust = 0),
          axis.text = element_text(colour = MUTED), axis.title = element_text(colour = INK),
          plot.background = element_rect(fill = "white", colour = NA))

pool[, share := n / sum(n)]
p1 <- ggplot(pool, aes(bin, share)) +
    geom_vline(xintercept = 0, colour = MUTED, linewidth = 0.4) +
    geom_col(fill = "#2a78d6", width = 0.1) +
    scale_y_continuous(labels = function(x) paste0(round(100 * x, 1), "%")) +
    coord_cartesian(xlim = c(-3, 3)) +
    labs(x = "log2FC = log2((purified + 1) / (bulk + 1))", y = "Share of transcript x Tumour values") + th

p2 <- ggplot(per_t, aes(pi, median_abs_log2fc)) +
    geom_point(aes(shape = split), colour = "#2a78d6", size = 2.2, alpha = 0.85) +
    scale_shape_manual(values = c(train = 16, test = 17), name = NULL) +
    labs(x = "DeMixT tumour proportion (RNA)", y = "Median |log2FC| per Tumour") + th +
    theme(legend.position = if (length(unique(per_t$split)) > 1) "bottom" else "none")

p3 <- ggplot(mk, aes(gene, log2fc)) +
    geom_hline(yintercept = 0, colour = MUTED, linewidth = 0.4) +
    geom_boxplot(fill = "#cde2fb", colour = "#2a78d6", outlier.size = 0.6, width = 0.7) +
    facet_grid(~ group, scales = "free_x", space = "free_x") +
    labs(x = NULL, y = "log2FC purified vs bulk (per Tumour)") + th +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = INK))

cairo_pdf(out_pdf, width = 10, height = 5, onefile = TRUE); print(p1); print(p2); print(p3); invisible(dev.off())
ggsave(sub("\\.pdf$", "_distribution.png", out_pdf), p1, width = 7, height = 4.5, dpi = 150)
ggsave(sub("\\.pdf$", "_per_tumour.png",   out_pdf), p2, width = 6, height = 4.5, dpi = 150)
ggsave(sub("\\.pdf$", "_markers.png",      out_pdf), p3, width = 10, height = 5, dpi = 150)
cat("->", out_pdf, "\n")
