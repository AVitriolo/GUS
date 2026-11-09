filterCov <- function(counts, minCov, minSamples=90){
	nb1 <- nrow(counts)
	to.keep <- which(apply(counts,1,minCov=minCov,minSamples=minSamples,function(x,minCov,minSamples){ sum(x>=minCov) })>=minSamples)
	message(paste(nb1-length(to.keep),"rows out of",nb1,"were removed [COVERAGE]."))
	gc()
	return(rownames(counts)[to.keep])
}