source("workflow/scripts/helpers/load_h5_rse.R")
source("workflow/scripts/helpers/process_beta_vals.R")
source("workflow/scripts/helpers/process_cov_vals.R")

library(future)
library(future.apply)

options(scipen=999)                                                                    #  unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                     #  read args

input_dir                                 <- args$input_dir                                                    #  resources/cpgea_wgbs_with_coverage_hg38/
input_dir_TMRs                           <- args$input_path_TMRs
input_path_K_closest                      <- args$input_path_K_closest
input_path_counts                         <- args$input_path_counts
sample_type                               <- args$sample_type
distance                                  <- as.numeric(args$distance)
min_CpG                                   <- as.numeric(args$min_CpG)
minCov                                    <- as.numeric(args$minCov)
leftCount_beta                            <- as.numeric(args$leftCount_beta)
rightCount_beta                           <- as.numeric(args$rightCount_beta)
minSamples_beta                           <- as.integer(args$minSamples_beta)
output_path_ff_selected_TxIDs             <- args$output_path_ff_selected_TxIDs
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

#### LOAD counts

counts <- read.table(input_path_counts, head=T)

#### Filter top K CpGs

distance.minCpGs.filt <- top_K_CpGs %>% 
                            dplyr::filter(.data$distance <= .env$distance) %>%                 # filter by distance .env here is needed because otherwise I compare the column with itself
                            dplyr::group_by(TxID) %>% 
                            dplyr::mutate(numCpGs = dplyr::n()) %>% 
                            dplyr::filter(numCpGs >= min_CpG)                                 # filter by numCpGs

TxIDs <- unique(distance.minCpGs.filt$TxID) # save

options(future.globals.maxSize = 4 * 1024^3)
plan(multicore, workers = 20)
future_lapply(X = TxIDs, FUN = function(id){                                                 # to parallelize
    
    counts_by_TxID <- as.data.frame(t(counts[id, grep(paste0("^",sample_type), colnames(counts))]))

    CpGs_by_TxID <- distance.minCpGs.filt %>% 
                        dplyr::filter(.data$TxID == id) %>%                                     # get CpGs per Tx
                        dplyr::pull(CpGID)

    idxs <- match(CpGs_by_TxID, GenomicRanges::mcols(CpGs)$CpGID)                               # retrieve the rowRanges indexes

    CpGs_by_TxID <- CpGs[idxs]                                                                       # subset CpGs
    
    cov_by_TxID <- as.data.frame(cov[idxs, grep(paste0("^",sample_type), colnames(cov))])            # subset CpGs and samples
    rownames(cov_by_TxID) <- GenomicRanges::mcols(CpGs_by_TxID)$CpGID

    beta_by_TxID <- as.data.frame(beta[idxs, grep(paste0("^",sample_type), colnames(beta))])         # subset CpGs and samples
    rownames(beta_by_TxID) <- GenomicRanges::mcols(CpGs_by_TxID)$CpGID

    CpGs_with_enough_coverage <- process_cov_vals(cov_by_TxID, minCov, minSamples_beta)              # get covered CpGs

    beta_by_TxID <- process_beta_vals(beta_by_TxID, leftCount_beta, rightCount_beta, minSamples_beta)     # filter by beta
    beta_by_TxID <- beta_by_TxID[rownames(beta_by_TxID) %in% CpGs_with_enough_coverage,]

    CpGs_by_TxID <- CpGs_by_TxID[which(GenomicRanges::mcols(CpGs_by_TxID)$CpGID %in% rownames(beta_by_TxID))]
    
    beta_by_TxID <- beta_by_TxID[match(GenomicRanges::mcols(CpGs_by_TxID)$CpGID,rownames(beta_by_TxID)),]
    
    coord_order <- GenomicRanges::mcols(GenomicRanges::sort(CpGs_by_TxID))$CpGID    
    
    gc()

    beta_by_TxID <- as.data.frame(t(beta_by_TxID))

    beta_by_TxID$sample <- rownames(beta_by_TxID)
    counts_by_TxID$sample <- rownames(counts_by_TxID)
    xgb_input <- merge(beta_by_TxID, counts_by_TxID, by = "sample")
    rownames(xgb_input) <- xgb_input$sample; xgb_input$sample <- NULL

    sd_per_value <- sapply(xgb_input[,1:(ncol(xgb_input)-1)], sd, na.rm = TRUE)
    most_variable_CpGs <- names(sd_per_value[order(sd_per_value, decreasing=T)][1:min(length(sd_per_value), (nrow(xgb_input)-1))])
    xgb_input <- xgb_input[,c(most_variable_CpGs, id)]

    xgb_input_wo_target <- xgb_input[,which(!(colnames(xgb_input) %in% id))]
    
    coord_order <- coord_order[coord_order %in% colnames(xgb_input_wo_target)]
    xgb_input_wo_target <- xgb_input_wo_target[, coord_order]

    correlations_table <- cor(xgb_input_wo_target, method = "spearman")
    
    basename <- strsplit(x = output_path_ff_selected_TxIDs, split = "/")[[1]][3]
    basename <- gsub(x = basename, pattern = "_ff_selected_TxIDs", replacement = "")
    
    output_path_corr <- paste0("data/corr/", basename, "_", id)
    output_path_xgb <- paste0("data/xgb/", basename, "_", id)

    write.table(correlations_table,
 		file=output_path_corr,
 		col.names = T,
 		row.names=T,
 		quote=F, 
 		sep="\t")

    write.table(xgb_input,
		file=output_path_xgb,
		col.names = T,
		row.names=T,
		quote=F, 
		sep="\t")

    return(id)

})
plan(sequential)

pdf(output_path_num_CpGs_plot)
hist(distance.minCpGs.filt$numCpGs, breaks = 15)
dev.off()

pdf(output_path_dist_filt_plot)
hist(log10(as.data.frame(distance.minCpGs.filt)$distance), breaks = 45)
dev.off()

writeLines(text = TxIDs, sep = "\n", con = output_path_ff_selected_TxIDs)

