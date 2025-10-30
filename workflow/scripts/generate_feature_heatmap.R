options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path                               <- args$input_path
output_path 	                         <- args$output_path

correlation_matrix <- read.table(input_path, row.names = 1, header = T)

ht <- ComplexHeatmap::Heatmap(
  as.matrix(correlation_matrix),
  name = "Correlation",
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", a[i, j]), x, y, gp = gpar(fontsize = 8, col = "black"))
  }
)

pdf(output_path, width = 6, height = 6)
ComplexHeatmap::draw(ht)
dev.off()
