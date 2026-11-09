suppressPackageStartupMessages({
  library(ggplot2)
})

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

puree_path <- args$puree_purity_m
coldata_path <- args$coldata
out_plot <- args$output_plot
out_table <- args$output_table

purity <- read.delim(puree_path)
colnames(purity)[1:2] <- c("sample", "purity")

cd <- readRDS(coldata_path)
cd <- data.frame(sample = rownames(cd), metastatis_site = cd$metastatis_site,
                 stringsAsFactors = FALSE)

merged <- merge(purity, cd, by = "sample")
merged <- merged[!is.na(merged$metastatis_site), ]

message(sprintf("%d of %d metastatic samples have a site annotation",
                nrow(merged), nrow(purity)))

# order sites by sample count, so the well-powered groups read first
site_n <- sort(table(merged$metastatis_site), decreasing = TRUE)
merged$metastatis_site <- factor(merged$metastatis_site, levels = names(site_n))
levels_with_n <- sprintf("%s\n(n=%d)", names(site_n), as.integer(site_n))
merged$site_label <- factor(levels_with_n[as.integer(merged$metastatis_site)],
                            levels = levels_with_n)

site_colors <- setNames(c("#2a78d6", "#eb6834", "#1baf7a", "#eda100")[seq_along(site_n)],
                        levels_with_n)

p <- ggplot(merged, aes(x = site_label, y = purity, fill = site_label)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.35, width = 0.6) +
  geom_jitter(aes(color = site_label), width = 0.15, height = 0, size = 2, alpha = 0.8) +
  scale_fill_manual(values = site_colors) +
  scale_color_manual(values = site_colors) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = NULL,
    y = "PUREE purity",
    title = "Metastatic sample purity by biopsy site",
    subtitle = sprintf("PUREE estimates, n = %d metastatic samples with site annotation", nrow(merged))
  ) +
  theme_bw() +
  theme(legend.position = "none", panel.grid.minor = element_blank())

ggsave(out_plot, p, width = 7, height = 5, dpi = 300)

summary_tbl <- aggregate(purity ~ metastatis_site, merged,
                        function(x) c(n = length(x), median = median(x),
                                      min = min(x), max = max(x)))
summary_tbl <- do.call(data.frame, summary_tbl)
print(summary_tbl)

write.table(merged[, c("sample", "purity", "metastatis_site")], out_table,
            sep = "\t", row.names = FALSE, quote = FALSE)
message(sprintf("Plot written to %s", out_plot))
