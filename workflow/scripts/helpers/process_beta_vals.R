source("workflow/scripts/helpers/filterBetas.R")

process_beta_vals <- function(beta_df, leftCount_beta = 0.25, rightCount_beta = 0.75, minSamples_beta){
    #rownames(beta_df) <- beta_df$CpG_ID
    #beta_df$CpG_ID <- NULL
    beta_df[is.na(beta_df)] <- 0
    beta_df.f <- filterBetas(beta_df,leftCount_beta=leftCount_beta, rightCount_beta=rightCount_beta, minSamples=minSamples_beta)
    gc()
    return(beta_df.f)
}