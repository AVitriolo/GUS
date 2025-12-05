source("workflow/scripts/helpers/load_h5_rse_mini.R")
source("workflow/scripts/helpers/process_beta_vals.R")
source("workflow/scripts/helpers/process_cov_vals.R")

options(scipen=999)                                                                                            #  unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                                             #  read args

input_dir_rse                             <- args$input_dir_rse                                                  #  resources/cpgea_wgbs_with_coverage_hg38/
input_dir_TMRs                            <- args$input_dir_TMRs
input_path_counts                         <- args$input_path_counts
TxID                                      <- args$TxID
sample_type                               <- args$sample_type
minCov                                    <- as.numeric(args$minCov)
leftCount_beta                            <- as.numeric(args$leftCount_beta)
rightCount_beta                           <- as.numeric(args$rightCount_beta)
minSamples_beta                           <- as.integer(args$minSamples_beta)
output_path_xgb                           <- args$output_path_xgb
output_path_corr                          <- args$output_path_corr

`%>%` <- magrittr::`%>%`

#### LOAD CpGs

input_dir_rse_TxID <- grep(pattern = TxID, list.files(input_dir_rse))

h5_list <- load_h5_rse_mini(input_dir_rse_TxID)
rse_by_TxID <- h5_list$rse
CpGs_by_TxID <- h5_list$CpGs
beta_by_TxID <- h5_list$beta
cov_by_TxID <- h5_list$cov

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

counts_by_TxID <- as.data.frame(t(counts[TxID, grep(paste0("^",sample_type), colnames(counts))]))


CpGs_with_enough_coverage <- process_cov_vals(cov_by_TxID, minCov, minSamples_beta)                        # get covered CpGs

beta_by_TxID <- process_beta_vals(beta_by_TxID, leftCount_beta, rightCount_beta, minSamples_beta)          # filter by beta
beta_by_TxID <- beta_by_TxID[rownames(beta_by_TxID) %in% CpGs_with_enough_coverage,]                       # filter by coverage
    
CpGs_by_TxID <- CpGs_by_TxID[which(GenomicRanges::mcols(CpGs_by_TxID)$CpGID %in% rownames(beta_by_TxID))]
    
beta_by_TxID <- beta_by_TxID[match(GenomicRanges::mcols(CpGs_by_TxID)$CpGID,rownames(beta_by_TxID)),]
    
coord_order <- GenomicRanges::mcols(GenomicRanges::sort(CpGs_by_TxID))$CpGID
    
gc()

TMRs_by_TxID <- TMRs.gr[GenomicRanges::mcols(TMRs.gr)$TxID == TxID]
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

non_CpGs <- colnames(xgb_input)[grep(pattern="^CpG", colnames(xgb_input), invert = TRUE)]                 # greps tmrID and TxID, they both follow ^id pattern

sd_per_value <- sapply(xgb_input[,!(colnames(xgb_input) %in% non_CpGs)], sd, na.rm = TRUE)
most_variable_CpGs <- names(sd_per_value[order(sd_per_value, decreasing=T)][1:min(length(sd_per_value), (nrow(xgb_input)-length(non_CpGs)))])
xgb_input <- xgb_input[,c(most_variable_CpGs, non_CpGs)]

xgb_input_wo_target <- xgb_input[,which(!(colnames(xgb_input) %in% TxID))]

coord_order <- coord_order[coord_order %in% colnames(xgb_input_wo_target)]
xgb_input_wo_target <- xgb_input_wo_target[, coord_order]

correlations_table <- cor(xgb_input_wo_target, method = "spearman")

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