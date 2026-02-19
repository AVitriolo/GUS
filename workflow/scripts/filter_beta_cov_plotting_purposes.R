source("workflow/scripts/helpers/load_h5_rse.R")
source("workflow/scripts/helpers/process_beta_vals.R")
source("workflow/scripts/helpers/process_cov_vals.R")

options(scipen=999)                                                                    #  unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                     #  read args

input_path_CpG_mat                        <- args$input_path_CpG_mat                       #  resources/cpgea_wgbs_with_coverage_hg38/
input_path_filtered_top_K                 <- args$input_path_filtered_top_K

sample_type                               <- args$sample_type
minCov                                    <- as.numeric(args$minCov)
leftCount_beta                            <- as.numeric(args$leftCount_beta)
rightCount_beta                           <- as.numeric(args$rightCount_beta)
minSamples_beta                           <- as.integer(args$minSamples_beta)

output_path                               <- args$output_path

`%>%` <- magrittr::`%>%`

#### LOAD CpGs

h5_list <- load_h5_rse(input_path_CpG_mat)
rse <- h5_list$rse
CpGs <- h5_list$CpGs
beta <- h5_list$beta
cov <- h5_list$cov

#### LOAD K closest

top_K_CpGs <- data.table::fread(input_path_filtered_top_K)
colnames(top_K_CpGs) <- c("TxID", "CpGID", "distance")

#### RETAIN only filtered topK CpGs

CpGs_to_keep_index <- which(GenomicRanges::mcols(CpGs)$CpGID %in% top_K_CpGs$CpGID)

CpGs.f <- CpGs[CpGs_to_keep_index]

beta.f <- beta[CpGs_to_keep_index, grep(paste0("^",sample_type), colnames(beta))]
rownames(beta.f) <- GenomicRanges::mcols(CpGs.f)$CpGID

print(colnames(beta.f))

cov.f  <- cov[CpGs_to_keep_index, grep(paste0("^",sample_type), colnames(cov))]
rownames(cov.f) <- GenomicRanges::mcols(CpGs.f)$CpGID

#### LOAD counts

CpGs_with_enough_coverage <- process_cov_vals(cov.f, minCov, minSamples_beta)              # get covered CpGs

beta.f <- process_beta_vals(beta.f, leftCount_beta, rightCount_beta, minSamples_beta)      # filter by beta
beta.f <- beta.f[rownames(beta.f) %in% CpGs_with_enough_coverage,, drop = FALSE]

top_K_CpGs.f <- top_K_CpGs[top_K_CpGs$CpGID %in% rownames(beta.f), c("CpGID", "TxID", "distance")]

write.table(x = top_K_CpGs.f, sep = "\t", file = output_path, col.names = T, quote = F, row.names = F)