suppressPackageStartupMessages({
  library(ggplot2)
})

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

pames_path <- args$pames_purity
out_plot <- args$output_plot

df <- read.delim(pames_path)
df$Condition <- factor(df$Condition, levels = c("Normal", "Tumour"))

medians <- aggregate(Purity ~ Condition, df, median)

group_colors <- c(Normal = "#2a78d6", Tumour = "#eb6834")

p <- ggplot(df, aes(x = Purity, fill = Condition, color = Condition)) +
  geom_density(alpha = 0.35, linewidth = 0.8) +
  geom_vline(data = medians, aes(xintercept = Purity, color = Condition),
            linetype = "dashed", linewidth = 0.6, show.legend = FALSE) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    x = "PAMES purity",
    y = "Density",
    title = "PAMES purity distribution by sample group",
    subtitle = sprintf("Medians: Normal=%.2f, Tumour=%.2f",
                       medians$Purity[medians$Condition == "Normal"],
                       medians$Purity[medians$Condition == "Tumour"]),
    fill = "Group", color = "Group"
  ) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

ggsave(out_plot, p, width = 7, height = 5, dpi = 300)
message(sprintf("Plot written to %s", out_plot))
