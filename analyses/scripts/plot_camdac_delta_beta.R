# =============================================================================
# plot_camdac_delta_beta.R
# How much CAMDAC changed the CPGEA Tumour betas: delta = m_t (purified) - m
# (bulk), over every CpG CAMDAC kept, all Tumours. Each pure.csv.gz is reduced
# to delta histogram counts + per-Tumour summaries (diploid CpGs vs CpGs in a
# called CNA peak), so the full tables are never held together.
#   page 1: pooled delta distribution, diploid vs CNA
#   page 2: per Tumour, median |delta| and share clipped to 0/1 vs purity
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(parallel) })

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
camdac  <- args$camdac_dir
n_cores <- as.integer(if (is.null(args$n_cores)) 4 else args$n_cores)
out_pdf <- args$out_pdf

files <- list.files(camdac, "\\.pure\\.csv\\.gz$", recursive = TRUE, full.names = TRUE)
BREAKS <- seq(-1, 1, by = 0.01)

one <- function(f) {
    p <- fread(f, select = c("m", "m_t", "m_t_raw", "CG_CN", "p"))
    # CNA = copy number differs from the sample baseline (2 in CPGEA, the ploidy in the metastases)
    p[, `:=`(delta = m_t - m, cna = abs(CG_CN - as.numeric(names(which.max(table(CG_CN))))) > 1e-6)]
    hist_of <- function(d) as.numeric(table(cut(d, BREAKS, include.lowest = TRUE)))
    list(
        hist = data.table(bin = head(BREAKS, -1) + 0.005,
                          diploid = hist_of(p[cna == FALSE, delta]),
                          cna     = hist_of(p[cna == TRUE,  delta])),
        stat = data.table(sample = sub("\\..*", "", basename(f)), purity = p$p[1],
                          n_cpgs = nrow(p), frac_cna = mean(p$cna),
                          median_absdelta = median(abs(p$delta)),
                          median_absdelta_diploid = median(abs(p[cna == FALSE, delta])),
                          median_absdelta_cna = if (any(p$cna)) median(abs(p[cna == TRUE, delta])) else NA_real_,
                          mean_delta = mean(p$delta),
                          frac_clipped = mean(p$m_t_raw < 0 | p$m_t_raw > 1)))
}

res   <- mclapply(files, one, mc.cores = n_cores)
stats <- rbindlist(lapply(res, `[[`, "stat"))
hist  <- rbindlist(lapply(res, `[[`, "hist"))[, .(diploid = sum(diploid), cna = sum(cna)), by = bin]
fwrite(stats, sub("\\.pdf$", "_per_tumour.tsv", out_pdf), sep = "\t")
fwrite(hist,  sub("\\.pdf$", "_histogram.tsv", out_pdf), sep = "\t")   # re-plot without re-reading

tot <- function(col) sum(hist[[col]])
cum_med <- function(col) { w <- hist[order(abs(bin))]; w$bin[which(cumsum(w[[col]]) >= tot(col) / 2)[1]] }
cat(sprintf("%d Tumours, %s CpG x Tumour values (%s in CNA peaks)\n", nrow(stats),
            format(tot("diploid") + tot("cna"), big.mark = ","), format(tot("cna"), big.mark = ",")))
cat(sprintf("median |delta|: diploid ~%.3f | CNA ~%.3f | mean delta across Tumours %.3f\n",
            abs(cum_med("diploid")), abs(cum_med("cna")), mean(stats$mean_delta)))
cat(sprintf("per Tumour: median |delta| %.3f (range %.3f-%.3f) | clipped %.1f%% (range %.1f-%.1f%%)\n",
            median(stats$median_absdelta), min(stats$median_absdelta), max(stats$median_absdelta),
            100 * median(stats$frac_clipped), 100 * min(stats$frac_clipped), 100 * max(stats$frac_clipped)))

long <- melt(hist, id.vars = "bin", variable.name = "group", value.name = "n")
long[, share := n / sum(n), by = group]
long[, group := factor(group, levels = c("diploid", "cna"),
                       labels = c(sprintf("Baseline copy number: %s values", format(tot("diploid"), big.mark = ",")),
                                  sprintf("CpGs in called CNA peaks: %s values", format(tot("cna"), big.mark = ","))))]

INK <- "#3d3d3a"; MUTED <- "#8a8984"; GRID <- "#e6e5df"
th <- theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
          strip.text = element_text(face = "bold", colour = INK, hjust = 0),
          axis.text = element_text(colour = MUTED), axis.title = element_text(colour = INK),
          plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = MUTED),
          plot.background = element_rect(fill = "white", colour = NA))

p1 <- ggplot(long, aes(bin, share)) +
    geom_vline(xintercept = 0, colour = MUTED, linewidth = 0.4) +
    geom_col(fill = "#2a78d6", width = 0.01) +
    facet_wrap(~ group, ncol = 2, scales = "free_y") +
    scale_y_continuous(labels = function(x) paste0(round(100 * x, 1), "%")) +
    coord_cartesian(xlim = c(-0.6, 0.6)) +
    labs(x = "Δβ = CAMDAC purified tumour beta (m_t) - bulk beta (m)", y = "Share of CpG x Tumour values") + th

sl <- rbind(data.table(purity = stats$purity, value = stats$median_absdelta, metric = "Median |Δβ| per Tumour"),
            data.table(purity = stats$purity, value = 100 * stats$frac_clipped, metric = "CpGs clipped to 0/1 (%)"))
p2 <- ggplot(sl, aes(purity, value)) +
    geom_point(colour = "#2a78d6", size = 2.2, alpha = 0.85) +
    facet_wrap(~ metric, ncol = 2, scales = "free_y") +
    labs(x = "PAMES purity", y = NULL,
         title = "Per Tumour: the lower the purity, the larger the correction",
         subtitle = sprintf("%d Tumours. Spearman with purity: median |Δβ| %.2f, clipped %.2f", nrow(stats),
                            cor(stats$purity, stats$median_absdelta, method = "s"),
                            cor(stats$purity, stats$frac_clipped, method = "s"))) + th

cairo_pdf(out_pdf, width = 11, height = 5, onefile = TRUE); print(p1); print(p2); invisible(dev.off())
ggsave(sub("\\.pdf$", "_distribution.png", out_pdf), p1, width = 11, height = 5, dpi = 150)
cat("->", out_pdf, "\n")
