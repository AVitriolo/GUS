filterCov <- function(counts, minCov, minSamples=90){
	nb1 <- nrow(counts)
	counts <- counts[which(apply(counts,1,minCov=minCov,minSamples=minSamples,function(x,minCov,minSamples){ sum(x>=minCov) })>=minSamples),] 
	message(paste(nb1-nrow(counts),"rows out of",nb1,"were removed [COVERAGE]."))
	return(counts)
}