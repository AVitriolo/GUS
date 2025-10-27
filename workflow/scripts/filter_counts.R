source("workflow/scripts/helpers/process_expr_vals.R")

options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_expr_values_path                   <- args$input_expr_values_path
minCount_expr							 <- as.numeric(args$minCount_expr)
minSamples_expr							 <- as.numeric(args$minSamples_expr)
output_path 	                         <- args$output_path

ExprValues.tmm <- process_expr_vals(input_expr_values_path, minCount_expr, minSamples_expr)

write.table(ExprValues.tmm,
		file=output_path,
		col.names = T,
		row.names=T,
		quote=F, 
		sep="\t")