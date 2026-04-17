options(scipen=999)                                                                    #  unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                     #  read args

input_path                                <- args$input_path                           
TxID                                      <- args$TxID
output_path                               <- args$output_path


xgb_input <- read.table(input_path, header = TRUE, sep = "\t")
expr <- xgb_input[,TxID]
cpgs <- xgb_input[,  !(colnames(xgb_input) %in% TxID), drop = FALSE]

cor_stat <- apply(X = cpgs, 
				  MARGIN = 2, 
				  FUN = function(cpg){
							test <- cor.test(cpg, expr)
							stat <- test$statistic
							}
						)

cor_pval <- apply(X = cpgs, 
				  MARGIN = 2, 
				  FUN = function(cpg){
							test <- cor.test(cpg, expr)
							pval <- test$p.value
							}
						)

names(cor_stat) == names(cor_pval)

cor_results_df <- data.frame(CpG = names(cor_stat), 
							 correlation = as.numeric(cor_stat), 
							 pvalue = as.numeric(cor_pval))

write.table(cor_results_df, file = output_path, sep = "\t", row.names = FALSE, quote = FALSE)

