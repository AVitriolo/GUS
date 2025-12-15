options(scipen=999)                                                                    #  unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                     #  read args

input_path                                <- args$input_path                                                #  resources/cpgea_wgbs_with_coverage_hg38/
TxID                                      <- args$TxID
output_path                               <- args$output_path


xgb_input <- read.table(input_path, header = TRUE, sep = "\t")
expr <- xgb_input[,TxID]
cpgs <- xgb_input[,  !(colnames(xgb_input) %in% TxID), drop = FALSE]
cor_results <- apply(cpgs, 2, function(cpg) cor(cpg, expr))
cor_results_df <- data.frame(CpG = names(cor_results), correlation = as.numeric(cor_results))
write.table(cor_results_df, file = output_path, 
	                sep = "\t", row.names = FALSE, quote = FALSE)

