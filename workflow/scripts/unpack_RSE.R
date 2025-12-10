source("workflow/scripts/helpers/load_h5_rse.R")

options(scipen=999)                                                                                            #  unable scientific notation           

BiocParallel::register(BiocParallel::MulticoreParam(workers = 8))
args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                                             #  read args

input_dir                                 <- args$input_dir                                                    #  resources/cpgea_wgbs_with_coverage_hg38/
input_path_filtered_CpGs                  <- args$input_path_filtered_CpGs
sample_type                               <- args$sample_type
output_path_rse                           <- args$output_path_rse

`%>%` <- magrittr::`%>%`

#### LOAD CpGs

h5_list <- load_h5_rse(input_dir)
rse <- h5_list$rse
CpGs <- h5_list$CpGs
beta <- h5_list$beta
cov <- h5_list$cov

#### LOAD Top K filtered CpGs

CpGs_by_TxID <- data.table::fread(input_path_filtered_CpGs)
colnames(CpGs_by_TxID) <- c("TxID", "CpGID")
CpGs_by_TxID <- CpGs_by_TxID %>% dplyr::pull(CpGID)

idxs <- match(CpGs_by_TxID, GenomicRanges::mcols(CpGs)$CpGID)                                              # retrieve the rowRanges indexes
CpGs_by_TxID <- CpGs[idxs]                                                                                 # subset CpGs

cov_by_TxID <- as.data.frame(cov[idxs, grep(paste0("^",sample_type), colnames(cov))])                      # subset CpGs and samples
rownames(cov_by_TxID) <- GenomicRanges::mcols(CpGs_by_TxID)$CpGID	

beta_by_TxID <- as.data.frame(beta[idxs, grep(paste0("^",sample_type), colnames(beta))])                   # subset CpGs and samples
rownames(beta_by_TxID) <- GenomicRanges::mcols(CpGs_by_TxID)$CpGID
    
gc()

rse_by_TxID <- SummarizedExperiment::SummarizedExperiment(assays = list(beta = beta_by_TxID, cov = cov_by_TxID), rowRanges = CpGs_by_TxID)
HDF5Array::saveHDF5SummarizedExperiment(x = rse_by_TxID, dir = output_path_rse)
