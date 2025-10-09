betaValues <- load...

rows = CpG
cols = Samples

1st : are the rownames == CpG_IDs? NO, extra column

ExprValues <- load ... # kallisto

rows = TxIDs
cols = Samples

library(edgeR)
donorm <- function (dataset, method = "TMM", returnCPM = FALSE)
{
    method <- match.arg(method, c("TMM", "RLE", "upperquartile",
        "quantile", "linear", "none"))
    if (method == "none")
        return(dataset)
    if (method == "linear") {
        nf <- getLinearNormalizers(dataset)
        return(t(t(dataset)/nf))
    }
    if (method == "quantile") {
        return(preprocessCore::normalize.quantiles(dataset))
    }
    d <- DGEList(counts = dataset, group = colnames(dataset))
    d <- calcNormFactors(d, method = method)
    e <- cpm(d, normalized.lib.sizes = T)
    if (!returnCPM)
        e <- e * mean(d$samples$lib.size)/1e+06
    return(as.data.frame(e))
}

# evaluate the option of log(donorm(EsprValues)

ExprValues.tmm <- donorm(ExprValues)

TransBeta <- t(betaValues)

CpG_Tx_distance <- load ...

names(CpG_Tx_distance) <- c("CpGID","TxID","distance")

for each TxIDs in CpG_Tx_distance {
	CpGs <- CpG_Tx_distance[which(CpG_Tx_distance$TxID == i),]
	Betas <- TransBeta[,which(colnames(TransBeta) %in% CpGs$CpGID)]
	samples <- intersect(rownames(Betas),colnames(ExprValues.tmm))
	Exprs <- t(ExprValues.tmm)
	Exprs <- Exprs[,i]
	Exprs <- as.data.frame(Exprs)
	Betas <- as.data.frame(Betas)
	colnames(Exprs) <- i
	Betas <- Betas[which(rownames(Betas) %in% samples),]
	Exprs <- Exprs[which(rownames(Exprs) %in% samples),]
	Betas$sample <- rownames(Betas)
	Exprs$sample <- rownames(Exprs)
	# ? rownames(Exprs) <- samples
	BetaExpr <- merge(Betas,Exprs,by="sample")
	rownames(BetaExpr) <- BetaExpr$sample
	BetaExpr$sample <- NULL
	BetaExpr[is.na(BetaExpr)] <- 0
	# Calcola Covarianza CpG_varz <- var(BetaExpr[,1:(ncol(BetaExpr)-1)])
	# for each column, lapply, var
done

