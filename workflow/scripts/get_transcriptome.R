source("workflow/scripts/helpers/process_expr_vals.R")

options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_filtered_expr_values_path          <- args$input_filtered_expr_values_path
TxID									 <- args$TxID
sample_type								 <- args$sample_type
output_path 	                         <- args$output_path

FilteredExprValues <- read.table(input_filtered_expr_values_path,head=T)

if(sample_type == "T"){

FilteredExprValues <- FilteredExprValues[, grep("^T", colnames(FilteredExprValues))]
FilteredExprValues <- as.data.frame(FilteredExprValues)

} else if (sample_type == "N"){

FilteredExprValues <- FilteredExprValues[, grep("^N", colnames(FilteredExprValues))]
FilteredExprValues <- as.data.frame(FilteredExprValues)

}

Exprs <- t(FilteredExprValues)
Exprs <- Exprs[,TxID]
Exprs <- as.data.frame(Exprs)
colnames(Exprs) <- TxID

write.table(Exprs,
		file=output_path,
		col.names = T,
		row.names=T,
		quote=F, 
		sep="\t")


