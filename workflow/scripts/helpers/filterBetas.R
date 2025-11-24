filterBetas <- function(counts, leftCount_beta=0.25, rightCount_beta=0.75, minSamples=90){
	nb1 <- nrow(counts)
	left_CpGs <- counts[which(apply(counts,1,leftCount_beta=leftCount_beta,minSamples=minSamples,function(x,leftCount_beta,minSamples){ sum(x<=leftCount_beta) })>=minSamples),] 
	right_CpGs <- counts[which(apply(counts,1,rightCount_beta=rightCount_beta,minSamples=minSamples,function(x,rightCount_beta,minSamples){ sum(x>=rightCount_beta) })>=minSamples),] 
	right_CpGs <- right_CpGs[-which(rownames(right_CpGs) %in% rownames(left_CpGs)),]
	counts <- rbind(left_CpGs,right_CpGs)
	rm(left_CpGs, right_CpGs)
	message(paste(nb1-nrow(counts),"rows out of",nb1,"were removed [BETAS]."))
	gc()
	return(counts)
}
