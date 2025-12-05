load_h5_rse_mini <- function(input_dir){
    
    rse <- HDF5Array::loadHDF5SummarizedExperiment(input_dir)
    
    CpGs <- rse@rowRanges
    beta <-SummarizedExperiment::assay(rse, "beta")
    cov <- SummarizedExperiment::assay(rse, "cov")

    out.list <- list(
        rse = rse,
        CpGs = CpGs,
        beta = beta,
        cov = cov
    )

    return(out.list)
    
}