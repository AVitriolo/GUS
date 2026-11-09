# =============================================================================
# select_demixt_grid.R
# Choose the DeMixT setting (nspikein x ngene) from the TRAIN grid, following
# the DeMixT vignette: the setting whose tumour proportions have the largest
# average Spearman correlation with those of every other setting. Also reports,
# per setting, the correlation with PAMES and how many Tumours sit at the
# upper plateau (pi >= 0.9) - for judgement, not for the selection itself.
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
grid_dir <- args$grid_dir
pames    <- args$pames
out_tsv  <- args$out_tsv
out_png  <- args$out_png

files <- list.files(grid_dir, "^demixt_pi_train_ns[0-9]+_ng[0-9]+\\.tsv$", full.names = TRUE)
pis   <- lapply(files, fread)
names(pis) <- sub("^demixt_pi_train_(.*)\\.tsv$", "\\1", basename(files))
stopifnot(all(sapply(pis, function(p) identical(p$sample, pis[[1]]$sample))))

M  <- sapply(pis, `[[`, "pi_tumour"); rownames(M) <- pis[[1]]$sample
pm <- fread(pames); setnames(pm, c("sample", "pames", "condition"))
pa <- pm$pames[match(rownames(M), pm$sample)]

C <- cor(M, method = "spearman")
res <- data.table(setting    = colnames(M),
                  nspikein   = as.integer(sub("ns([0-9]+)_ng.*", "\\1", colnames(M))),
                  ngene      = as.integer(sub(".*_ng([0-9]+)", "\\1", colnames(M))),
                  mean_cor_other_settings = (colSums(C) - 1) / (ncol(C) - 1),
                  spearman_pames = apply(M, 2, cor, y = pa, method = "spearman"),
                  median_pi  = apply(M, 2, median),
                  frac_pi_ge_0.9 = colMeans(M >= 0.9))
res <- res[order(-mean_cor_other_settings)]
res[, selected := seq_len(.N) == 1]
fwrite(res, out_tsv, sep = "\t")
print(res[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)], row.names = FALSE)
cat(sprintf("selected: %s\n", res$setting[1]))

INK <- "#3d3d3a"; MUTED <- "#8a8984"; GRID <- "#e6e5df"
res[, ngene_f := factor(ngene)]; res[, nspikein_f := factor(nspikein)]
long <- melt(res, id.vars = c("nspikein_f", "ngene_f", "selected"),
             measure.vars = c("mean_cor_other_settings", "spearman_pames", "frac_pi_ge_0.9"))
long[, variable := factor(variable, levels = c("mean_cor_other_settings", "spearman_pames", "frac_pi_ge_0.9"),
                          labels = c("Mean Spearman with other settings (selection rule)",
                                     "Spearman with PAMES", "Share of Tumours with pi >= 0.9"))]
p <- ggplot(long, aes(ngene_f, nspikein_f, fill = value)) +
    geom_tile(colour = "white", linewidth = 1) +
    geom_text(aes(label = sprintf("%.2f", value), fontface = ifelse(selected, "bold", "plain")),
              colour = INK, size = 3.3) +
    scale_fill_gradient(low = "#cde2fb", high = "#2a78d6", guide = "none") +
    facet_wrap(~ variable, ncol = 3) +
    labs(x = "Genes used for pi (ngene)", y = "Spike-ins (nspikein)",
         title = "DeMixT parameter grid on train Tumours",
         subtitle = sprintf("Selected (bold): %s. Each panel is shaded on its own scale", res$setting[1])) +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(), strip.text = element_text(face = "bold", colour = INK, hjust = 0),
          axis.text = element_text(colour = MUTED), axis.title = element_text(colour = INK),
          plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = MUTED),
          plot.background = element_rect(fill = "white", colour = NA))
ggsave(out_png, p, width = 13, height = 4.6, dpi = 150)
cat("->", out_png, "\n")
