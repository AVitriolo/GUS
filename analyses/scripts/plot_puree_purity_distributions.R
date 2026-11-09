suppressPackageStartupMessages({
  library(ggplot2)
})

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

n_path <- args$puree_n
t_path <- args$puree_t
m_path <- args$puree_m
out_plot <- args$output_plot

read_group <- function(path, label) {
  df <- read.delim(path)
  colnames(df)[1:2] <- c("sample", "purity")
  df$group <- label
  df[, c("sample", "purity", "group")]
}

all_purity <- rbind(
  read_group(n_path, "Normal"),
  read_group(t_path, "Tumour"),
  read_group(m_path, "Metastatic")
)
all_purity$group <- factor(all_purity$group, levels = c("Normal", "Tumour", "Metastatic"))

medians <- aggregate(purity ~ group, all_purity, median)

group_colors <- c(Normal = "#2a78d6", Tumour = "#eb6834", Metastatic = "#1baf7a")

p <- ggplot(all_purity, aes(x = purity, fill = group, color = group)) +
  geom_density(alpha = 0.35, linewidth = 0.8) +
  geom_vline(data = medians, aes(xintercept = purity, color = group),
            linetype = "dashed", linewidth = 0.6, show.legend = FALSE) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    x = "PUREE purity",
    y = "Density",
    title = "PUREE purity distribution by sample group",
    subtitle = sprintf("Medians: Normal=%.2f, Tumour=%.2f, Metastatic=%.2f",
                       medians$purity[medians$group == "Normal"],
                       medians$purity[medians$group == "Tumour"],
                       medians$purity[medians$group == "Metastatic"]),
    fill = "Group", color = "Group"
  ) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

ggsave(out_plot, p, width = 7, height = 5, dpi = 300)
message(sprintf("Plot written to %s", out_plot))
