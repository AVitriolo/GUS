compute_corr <- function(input_df){
    cor_mat <- cor(input_df)
    mean_abs_corr <- mean(abs(cor_mat))
    return(mean_abs_corr)
}