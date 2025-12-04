source("workflow/scripts/helpers/load_h5_rse.R")
source("workflow/scripts/helpers/process_beta_vals.R")
source("workflow/scripts/helpers/process_cov_vals.R")
#source("workflow/scripts/helpers/summarizeRegionMethylation.R")
#source("workflow/scripts/helpers/summarize_chunk_methylation.R")
#source("workflow/scripts/helpers/chunk_regions.R")

#library(future)
#library(MatrixGenerics)
#library(GenomeInfoDb)

options(scipen=999)                                                                                            #  unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                                             #  read args

input_dir                                 <- args$input_dir                                                    #  resources/cpgea_wgbs_with_coverage_hg38/
input_dir_TMRs                            <- args$input_dir_TMRs
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

##### LOAD TMRs

input_path_TMRs <- grep(pattern = paste0("^",sample_type), x = list.files(input_dir_TMRs), value = T)
input_path_TMRs <- paste0(input_dir_TMRs, "/", input_path_TMRs)

TMRs <- data.table::fread(input_path_TMRs, header = F)
colnames(TMRs) <- c("chr", "start", "end", "tmrID", "unk1", "unk2")
TMRs$unk1 <- NULL; TMRs$unk2 <- NULL 

TMRs.gr <- GenomicRanges::GRanges(TMRs)  #: input 3 di summarizeRegionMethylation "genomic_regions"
GenomicRanges::mcols(TMRs.gr)$TxID <- lapply(X = TMRs$tmrID, FUN = function(x){strsplit(x, split = "_")[[1]][1]})

#### Filter top K CpGs

distance.minCpGs.filt <- top_K_CpGs %>% 
                            dplyr::filter(.data$distance <= .env$distance) %>%                                 # filter by distance .env here is needed because otherwise I compare the column with itself
                            dplyr::group_by(TxID) %>% 
                            dplyr::mutate(numCpGs = dplyr::n()) %>% 
                            dplyr::filter(numCpGs >= min_CpG)                                                  # filter by numCpGs

TxIDs <- unique(distance.minCpGs.filt$TxID)


lapply(X = TxIDs, FUN = function(id){                                                                          # to parallelize
    
    counts_by_TxID <- as.data.frame(t(counts[id, grep(paste0("^",sample_type), colnames(counts))]))

    CpGs_by_TxID <- distance.minCpGs.filt %>% 
                        dplyr::filter(.data$TxID == id) %>%                                                    # get CpGs per Tx
                        dplyr::pull(CpGID)

    idxs <- match(CpGs_by_TxID, GenomicRanges::mcols(CpGs)$CpGID)                                              # retrieve the rowRanges indexes

    CpGs_by_TxID <- CpGs[idxs]                                                                                 # subset CpGs
    
    cov_by_TxID <- as.data.frame(cov[idxs, grep(paste0("^",sample_type), colnames(cov))])                      # subset CpGs and samples
    rownames(cov_by_TxID) <- GenomicRanges::mcols(CpGs_by_TxID)$CpGID

    beta_by_TxID <- as.data.frame(beta[idxs, grep(paste0("^",sample_type), colnames(beta))])                   # subset CpGs and samples
    rownames(beta_by_TxID) <- GenomicRanges::mcols(CpGs_by_TxID)$CpGID

    CpGs_with_enough_coverage <- process_cov_vals(cov_by_TxID, minCov, minSamples_beta)                        # get covered CpGs

    beta_by_TxID <- process_beta_vals(beta_by_TxID, leftCount_beta, rightCount_beta, minSamples_beta)          # filter by beta
    beta_by_TxID <- beta_by_TxID[rownames(beta_by_TxID) %in% CpGs_with_enough_coverage,]                       # filter by coverage
    
    CpGs_by_TxID <- CpGs_by_TxID[which(GenomicRanges::mcols(CpGs_by_TxID)$CpGID %in% rownames(beta_by_TxID))]
    
    beta_by_TxID <- beta_by_TxID[match(GenomicRanges::mcols(CpGs_by_TxID)$CpGID,rownames(beta_by_TxID)),]
    
    coord_order <- GenomicRanges::mcols(GenomicRanges::sort(CpGs_by_TxID))$CpGID
    
    gc()

    #TMRs_by_TxID <- TMRs.gr[GenomicRanges::mcols(TMRs.gr)$TxID == id]
    #TMRs_by_TxID_iCpGs <- IRanges::subsetByOverlaps(TMRs_by_TxID, CpGs_by_TxID, invert = FALSE)
    #iCpGs <- IRanges::subsetByOverlaps(CpGs_by_TxID, TMRs_by_TxID_iCpGs, invert = FALSE)
    #beta_by_TxID_iCpGs <- beta_by_TxID[rownames(beta_by_TxID) %in% GenomicRanges::mcols(iCpGs)$CpGID, ]
    #rse_iCpGs <-  SummarizedExperiment::SummarizedExperiment(assays = list(beta = beta_by_TxID_iCpGs), rowRanges = iCpGs)
    #avg_beta_by_TxID <- summarizeRegionMethylation(rse_iCpGs, 1, TMRs_by_TxID_iCpGs, genomic_region_names = GenomicRanges::mcols(TMRs_by_TxID_iCpGs)$tmrID)

    TMRs_by_TxID <- TMRs.gr[GenomicRanges::mcols(TMRs.gr)$TxID == id]
    TMRs_by_TxID_with_at_least_one_CpG <- IRanges::subsetByOverlaps(TMRs_by_TxID, CpGs_by_TxID, invert = FALSE)

    if(length(TMRs_by_TxID_with_at_least_one_CpG) > 0){

      beta_by_TxID_iCpGs <- list()

      for(idx in seq(length(TMRs_by_TxID_with_at_least_one_CpG))){
        
        TMR <- TMRs_by_TxID_with_at_least_one_CpG[idx]
        tmrID <- GenomicRanges::mcols(TMR)$tmrID
        iCpGs <- IRanges::subsetByOverlaps(CpGs_by_TxID, TMR, invert = FALSE)
        beta_by_TxID_by_TMR <- beta_by_TxID[rownames(beta_by_TxID) %in% GenomicRanges::mcols(iCpGs)$CpGID, ]
        avg_beta_by_TxID_by_TMR <- apply(X = beta_by_TxID_by_TMR, MARGIN = 2, FUN = mean)

        beta_by_TxID_iCpGs[[tmrID]] <- avg_beta_by_TxID_by_TMR
      }

      beta_by_TxID_iCpGs <- do.call(rbind, beta_by_TxID_iCpGs)
      
      oCpGs <- IRanges::subsetByOverlaps(CpGs_by_TxID, TMRs_by_TxID_with_at_least_one_CpG, invert = TRUE)
      beta_by_TxID_oCpGs <- beta_by_TxID[rownames(beta_by_TxID) %in% GenomicRanges::mcols(oCpGs)$CpGID, ]

      beta_by_TxID <- rbind(beta_by_TxID_oCpGs, beta_by_TxID_iCpGs)

      GenomicRanges::mcols(oCpGs)$ID <- GenomicRanges::mcols(oCpGs)$CpGID
      GenomicRanges::mcols(TMRs_by_TxID_with_at_least_one_CpG)$ID <- GenomicRanges::mcols(TMRs_by_TxID_with_at_least_one_CpG)$tmrID
      GenomicRanges::mcols(TMRs_by_TxID_with_at_least_one_CpG)$tmrID <- NULL
      GenomicRanges::mcols(oCpGs)$CpGID <- NULL
      GenomicRanges::mcols(TMRs_by_TxID_with_at_least_one_CpG)$TxID <- NULL

      CpGs_and_TMRs <- GenomicRanges::sort(c(oCpGs, TMRs_by_TxID_with_at_least_one_CpG))

      coord_order <- GenomicRanges::mcols(CpGs_and_TMRs)$ID

    } 

    beta_by_TxID <- as.data.frame(t(beta_by_TxID))
    beta_by_TxID$sample <- rownames(beta_by_TxID)

    counts_by_TxID$sample <- rownames(counts_by_TxID)
    xgb_input <- merge(beta_by_TxID, counts_by_TxID, by = "sample")
    rownames(xgb_input) <- xgb_input$sample; xgb_input$sample <- NULL

    non_CpGs <- colnames(xgb_input)[grep(pattern="^CpG", colnames(xgb_input), invert = TRUE)]                  # greps tmrID and TxID, they both follow ^id pattern
  
    sd_per_value <- sapply(xgb_input[,!(colnames(xgb_input) %in% non_CpGs)], sd, na.rm = TRUE)
    most_variable_CpGs <- names(sd_per_value[order(sd_per_value, decreasing=T)][1:(nrow(xgb_input)-length(non_CpGs))])
    xgb_input <- xgb_input[,c(most_variable_CpGs, non_CpGs)]

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

pdf(output_path_num_CpGs_plot)
hist(distance.minCpGs.filt$numCpGs, breaks = 15)
dev.off()

pdf(output_path_dist_filt_plot)
hist(log10(as.data.frame(distance.minCpGs.filt)$distance), breaks = 45)
dev.off()

writeLines(text = TxIDs, sep = "\n", con = output_path_ff_selected_TxIDs)

