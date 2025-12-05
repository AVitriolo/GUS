source("workflow/scripts/helpers/load_h5_rse.R")

options(scipen=999)                                                                                            #  unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                                             #  read args

input_dir                                 <- args$input_dir                                                    #  resources/cpgea_wgbs_with_coverage_hg38/
input_path_K_closest                      <- args$input_path_K_closest
input_path_selected_TxIDs                 <- args$input_path_selected_TxIDs
sample_type                               <- args$sample_type
distance                                  <- as.numeric(args$distance)
min_CpG                                   <- as.numeric(args$min_CpG)
output_dir                                <- args$output_dir
output_path_num_CpGs_plot                 <- args$output_path_num_CpGs_plot
output_path_dist_filt_plot                <- args$output_path_dist_filt_plot

`%>%` <- magrittr::`%>%`

#### LOAD CpGs

h5_list <- load_h5_rse(input_dir)
rse <- h5_list$rse
CpGs <- h5_list$CpGs
beta <- h5_list$beta
cov <- h5_list$cov

#### LOAD Top K CpGs

top_K_CpGs <- data.table::fread(input_path_K_closest)
colnames(top_K_CpGs) <- c("chr", "start", "end", "CpGID", "TxID", "distance")

#### LOAD TxIDs

TxIDs <- readLines(input_path_selected_TxIDs)

#### Filter top K CpGs

distance.minCpGs.filt <- top_K_CpGs %>% 
                            dplyr::filter(.data$distance <= .env$distance) %>%                                 # filter by distance .env here is needed because otherwise I compare the column with itself
                            dplyr::group_by(TxID) %>% 
                            dplyr::mutate(numCpGs = dplyr::n()) %>% 
                            dplyr::filter(numCpGs >= min_CpG)                                                  # filter by numCpGs

pdf(output_path_num_CpGs_plot)
hist(distance.minCpGs.filt$numCpGs, breaks = 15)
dev.off()

pdf(output_path_dist_filt_plot)
hist(log10(as.data.frame(distance.minCpGs.filt)$distance), breaks = 45)
dev.off()

lapply(X = TxIDs, FUN = function(id){                                                                          # to parallelize
    
    CpGs_by_TxID <- distance.minCpGs.filt %>% 
                        dplyr::filter(.data$TxID == id) %>%                                                    # get CpGs per Tx
                        dplyr::pull(CpGID)

    idxs <- match(CpGs_by_TxID, GenomicRanges::mcols(CpGs)$CpGID)                                              # retrieve the rowRanges indexes

    CpGs_by_TxID <- CpGs[idxs]                                                                                 # subset CpGs
    
    cov_by_TxID <- as.data.frame(cov[idxs, grep(paste0("^",sample_type), colnames(cov))])                      # subset CpGs and samples
    rownames(cov_by_TxID) <- GenomicRanges::mcols(CpGs_by_TxID)$CpGID

    beta_by_TxID <- as.data.frame(beta[idxs, grep(paste0("^",sample_type), colnames(beta))])                   # subset CpGs and samples
    rownames(beta_by_TxID) <- GenomicRanges::mcols(CpGs_by_TxID)$CpGID
    
    gc()

    rse_by_TxID <- SummarizedExperiment::SummarizedExperiment(assays = list(beta = beta_by_TxID, cov = cov_by_TxID), rowRanges = CpGs_by_TxID)
    HDF5Array::saveHDF5SummarizedExperiment(x = rse_by_TxID, dir = paste0(output_dir, "/", id))
    return(id)

})

