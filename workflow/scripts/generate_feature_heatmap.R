options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path                               <- args$input_path
output_path 	                           <- args$output_path

correlation_matrix <- read.table(input_path, row.names = 1, header = T)

mat <- as.matrix(correlation_matrix)

idxs <- order(unlist(lapply(X = rownames(mat), 
                            FUN = function(x){
                              as.numeric(strsplit(x, split = "_")[[1]][[2]])
                              })), 
                              decreasing = F)
mat <- mat[idxs, idxs]

col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

ht <- ComplexHeatmap::Heatmap(
  mat,
  name = "Correlation",
  cluster_rows = F,
  cluster_columns = F,
  column_names_gp = grid::gpar(fontsize = 4),
  row_names_gp = grid::gpar(fontsize = 4),
  col = col_fun
)

pdf(output_path, width = 6, height = 6)
ComplexHeatmap::draw(ht)
dev.off()
