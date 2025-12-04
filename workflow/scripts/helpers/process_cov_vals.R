source("workflow/scripts/helpers/filterCov.R")

process_cov_vals <- function(cov_df, minCov = 10, minSamples_beta){
    #rownames(cov_df) <- cov_df$CpG_ID 
    #cov_df$CpG_ID <- NULL
    cov_df[is.na(cov_df)] <- 0
    to.return <- filterCov(cov_df, minCov = minCov, minSamples=minSamples_beta)
    gc()
    return(to.return)
}