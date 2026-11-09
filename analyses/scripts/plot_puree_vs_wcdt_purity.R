suppressPackageStartupMessages({
  library(ggplot2)
})

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

puree_path <- args$puree_purity      # resources/purity/puree_purity_M.tsv
wcdt_path  <- args$wcdt_rdata        # ~/Lundberg_et_al/2022_11_04_WCDT_poppy.Rdata
out_plot   <- args$output_plot
out_table  <- args$output_table

e <- new.env()
load(wcdt_path, envir = e)
so <- get("SO_WCDT_WGS", envir = e)

wcdt <- data.frame(
  sample_id = sapply(so, function(x) x$sample_id),
  purity_wgs = sapply(so, function(x) if (is.null(x$purity)) NA else x$purity)
)
wcdt$base  <- sub("-PRO[0-9]*$|-BL$", "", wcdt$sample_id)
wcdt$is_bl <- grepl("-BL$", wcdt$sample_id)
# When a patient has multiple WGS timepoints, prefer the baseline (-BL) biopsy,
# since the RNA-seq/WGBS metastatic samples used elsewhere in this pipeline are
# single-timepoint and most closely match the baseline collection.
wcdt <- wcdt[order(wcdt$base, -wcdt$is_bl), ]
wcdt <- wcdt[!duplicated(wcdt$base), ]

puree <- read.delim(puree_path)
colnames(puree)[1] <- "sample_id"
colnames(puree)[2] <- "purity_puree"
puree$base <- sub("_", "-", puree$sample_id)

merged <- merge(puree[, c("base", "purity_puree")], wcdt[, c("base", "purity_wgs")], by = "base")

assertthat::assert_that(nrow(merged) > 0, msg = "No overlapping samples between PUREE and WCDT purity data.")

pearson  <- cor(merged$purity_puree, merged$purity_wgs, method = "pearson")
spearman <- cor(merged$purity_puree, merged$purity_wgs, method = "spearman")

message(sprintf("n = %d matched samples; Pearson r = %.3f; Spearman rho = %.3f",
                nrow(merged), pearson, spearman))

p <- ggplot(merged, aes(x = purity_wgs, y = purity_puree)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "steelblue") +
  coord_cartesian(xlim = c(0, NA), ylim = c(0, NA)) +
  labs(
    x = "WGS purity (PURPLE, Lundberg et al. WCDT)",
    y = "RNA-seq purity (PUREE)",
    title = "Metastatic (DTB) sample purity: PUREE vs WGS-PURPLE",
    subtitle = sprintf("n = %d, Pearson r = %.2f, Spearman rho = %.2f", nrow(merged), pearson, spearman)
  ) +
  theme_bw()

ggsave(out_plot, p, width = 6, height = 6, dpi = 300)
write.table(merged, out_table, sep = "\t", row.names = FALSE, quote = FALSE)

message(sprintf("Plot written to %s, table written to %s", out_plot, out_table))
