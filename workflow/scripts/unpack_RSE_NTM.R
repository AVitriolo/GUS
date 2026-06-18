source("workflow/scripts/helpers/load_h5_rse_integrated.R")

options(scipen=999)

BiocParallel::register(BiocParallel::MulticoreParam(workers = 8))

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

input_dir                <- args$input_dir
input_path_filtered_CpGs <- args$input_path_filtered_CpGs
output_path_rse          <- args$output_path_rse

`%>%` <- magrittr::`%>%`

# ---- load integrated SE ----
h5_list <- load_h5_rse_integrated(input_dir)
rse  <- h5_list$rse
CpGs <- h5_list$CpGs
M    <- h5_list$M

# ---- load per-TxID CpG list ----
CpGs_by_TxID <- data.table::fread(input_path_filtered_CpGs, header = FALSE)
colnames(CpGs_by_TxID) <- c("TxID", "CpGID")
CpGs_by_TxID <- CpGs_by_TxID %>% dplyr::pull(CpGID)

# ---- subset to TxID CpGs ----
idxs        <- match(CpGs_by_TxID, GenomicRanges::mcols(CpGs)$CpGID)
CpGs_subset <- CpGs[idxs]

M_by_TxID <- as.matrix(M[idxs, , drop = FALSE])

gc()

rse_by_TxID <- SummarizedExperiment::SummarizedExperiment(
    assays    = list(M = M_by_TxID),
    rowRanges = CpGs_subset,
    colData   = SummarizedExperiment::colData(rse)
)
rownames(rse_by_TxID) <- GenomicRanges::mcols(CpGs_subset)$CpGID

HDF5Array::saveHDF5SummarizedExperiment(
    rse_by_TxID,
    dir     = output_path_rse,
    replace = TRUE
)
