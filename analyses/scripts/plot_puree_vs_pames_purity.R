suppressPackageStartupMessages({
  library(ggplot2)
})

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

puree_path <- args$puree_purity   # resources/purity/puree_purity_T.tsv
pames_path <- args$pames_purity   # resources/purity/pames_purity_NT.tsv
out_plot   <- args$output_plot
out_table  <- args$output_table

puree <- read.delim(puree_path)
colnames(puree)[1] <- "sample_id"
colnames(puree)[2] <- "purity_puree"

pames <- read.delim(pames_path)
pames <- pames[pames$Condition == "Tumour", c("Sample", "Purity")]
colnames(pames) <- c("sample_id", "purity_pames")

merged <- merge(puree, pames, by = "sample_id")
assertthat::assert_that(nrow(merged) > 0, msg = "No overlapping tumour samples between PUREE and PAMES.")

pearson  <- cor(merged$purity_puree, merged$purity_pames, method = "pearson")
spearman <- cor(merged$purity_puree, merged$purity_pames, method = "spearman")

message(sprintf("n = %d matched tumour samples; Pearson r = %.3f; Spearman rho = %.3f",
                nrow(merged), pearson, spearman))

p <- ggplot(merged, aes(x = purity_pames, y = purity_puree)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "steelblue") +
  coord_cartesian(xlim = c(0, NA), ylim = c(0, NA)) +
  labs(
    x = "WGBS purity (PAMES)",
    y = "RNA-seq purity (PUREE)",
    title = "Primary tumor sample purity: PUREE vs PAMES",
    subtitle = sprintf("n = %d, Pearson r = %.2f, Spearman rho = %.2f", nrow(merged), pearson, spearman)
  ) +
  theme_bw()

ggsave(out_plot, p, width = 6, height = 6, dpi = 300)
write.table(merged, out_table, sep = "\t", row.names = FALSE, quote = FALSE)

message(sprintf("Plot written to %s, table written to %s", out_plot, out_table))
