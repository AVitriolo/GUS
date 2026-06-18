load_h5_rse_integrated <- function(input_h5,
                                    out_rowranges = gsub("\\.h5$", "_rowRanges.rds", input_h5),
                                    out_coldata   = gsub("\\.h5$", "_colData.rds",   input_h5)) {

    M        <- HDF5Array::HDF5Array(input_h5, "M")
    col_data <- readRDS(out_coldata)
    row_gr   <- readRDS(out_rowranges)

    rse <- SummarizedExperiment::SummarizedExperiment(
        assays    = list(M = M),
        rowRanges = row_gr,
        colData   = S4Vectors::DataFrame(col_data)
    )

    CpGs <- SummarizedExperiment::rowRanges(rse)
    GenomicRanges::mcols(CpGs)$CpGID <- paste0(
        GenomicRanges::seqnames(CpGs), "_",
        GenomicRanges::start(CpGs)
    )

    out.list <- list(
        rse  = rse,
        CpGs = CpGs,
        M    = M
    )

    return(out.list)
}
