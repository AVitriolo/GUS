load_h5_rse <- function(input_dir){
    
    rse <- HDF5Array::loadHDF5SummarizedExperiment(input_dir)
    
    CpGs <- rse@rowRanges
    #GenomicRanges::mcols(CpGs)$CpGID <- paste0("CpG_", 1:length(CpGs))
    GenomicRanges::mcols(CpGs)$CpGID <- paste0(GenomicRanges::seqnames(CpGs), "_", GenomicRanges::start(CpGs))

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
