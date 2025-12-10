source("workflow/scripts/helpers/process_expr_vals.R")

options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("edgeR")

options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_expr_values_path                   <- args$input_expr_values_path
minCount_expr							 <- as.numeric(args$minCount_expr)
minSamples_expr							 <- as.numeric(args$minSamples_expr)
output_path_counts 	                     <- args$output_path_counts
output_path_TxIDs						 <- args$output_path_TxIDs


ExprValues.tmm <- process_expr_vals(input_expr_values_path, minCount_expr, minSamples_expr)

write.table(ExprValues.tmm,
		file=output_path_counts,
		col.names = T,
		row.names=T,
		quote=F, 
		sep="\t")

writeLines(text = rownames(ExprValues.tmm), con = output_path_TxIDs, sep = "\n")
