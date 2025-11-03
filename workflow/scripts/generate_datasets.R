options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_transcriptome                      <- args$input_transcriptome
input_subset_beta                        <- args$input_subset_beta
TxID									 <- args$TxID
K                                        <- as.numeric(args$K)
output_path_xgb 	                     <- args$output_path_xgb
output_path_corr_table				 	 <- args$output_path_corr_table

Exprs <- read.table(input_transcriptome, head=T)
Betas <- t(read.table(input_subset_beta,head=T))

if(K >= nrow(Betas)){
	K <- nrow(Betas) - 1 # add mapping to closest value
} else {
	K <- K
}

samples <- intersect(rownames(Betas),rownames(Exprs))

Exprs <- Exprs[which(rownames(Exprs) %in% samples), , drop = F]
Betas <- as.data.frame(Betas[which(rownames(Betas) %in% samples),])

Betas$sample <- rownames(Betas)
Exprs$sample <- rownames(Exprs)

BetaExpr <- merge(Betas,Exprs,by="sample")

rownames(BetaExpr) <- BetaExpr$sample
BetaExpr$sample <- NULL

sd_per_value <- sapply(BetaExpr[,1:(ncol(BetaExpr)-1)], sd, na.rm = TRUE)
topKCpGs <- names(sd_per_value[order(sd_per_value, decreasing=T)][1:K])

correlations_table <- cor(BetaExpr[,which(colnames(BetaExpr) %in% topKCpGs)], method = "spearman")

Cols2keep <- cbind(topKCpGs,TxID)
BetaExpr.f <- BetaExpr[,which(colnames(BetaExpr) %in% Cols2keep)]

write.table(correlations_table,
 		file=output_path_corr_table,
 		col.names = T,
 		row.names=T,
 		quote=F, 
 		sep="\t")

write.table(BetaExpr.f,
		file=output_path_xgb,
		col.names = T,
		row.names=T,
		quote=F, 
		sep="\t")