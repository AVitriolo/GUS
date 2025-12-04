source("workflow/scripts/helpers/process_beta_vals.R")
source("workflow/scripts/helpers/process_cov_vals.R")

options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path                    <- args$input_path
output_path_beta              <- args$output_path_beta
output_path_CpGs              <- args$output_path_CpGs
sample_type                   <- args$sample_type
minCov                        <- as.numeric(args$minCov)
leftCount_beta                <- as.numeric(args$leftCount_beta)
rightCount_beta               <- as.numeric(args$rightCount_beta)
minSamples_beta               <- as.integer(args$minSamples_beta)

rse <- HDF5Array::loadHDF5SummarizedExperiment(input_path)

coord  <- SummarizedExperiment::rowRanges(rse)
beta_df <-SummarizedExperiment::assay(rse, "beta")
cov_df <- SummarizedExperiment::assay(rse, "cov")

rm(rse); gc()

if(sample_type == "T"){

beta_df <- beta_df[, grep("^T", colnames(beta_df))]; cov_df <- cov_df[, grep("^T", colnames(cov_df))]
beta_df <- as.data.frame(beta_df); cov_df <- as.data.frame(cov_df)
beta_df$CpG_ID = paste0("CpG_", 1:nrow(beta_df)); cov_df$CpG_ID = paste0("CpG_", 1:nrow(cov_df))

} else if (sample_type == "N"){

beta_df <- beta_df[, grep("^N", colnames(beta_df))]; cov_df <- cov_df[, grep("^N", colnames(cov_df))]
beta_df <- as.data.frame(beta_df); cov_df <- as.data.frame(cov_df)
beta_df$CpG_ID = paste0("CpG_", 1:nrow(beta_df)); cov_df$CpG_ID = paste0("CpG_", 1:nrow(cov_df))

}

CpGs_with_enough_coverage <- process_cov_vals(cov_df, minCov, minSamples_beta)
rm(cov_df); gc()
beta_df <- process_beta_vals(beta_df, leftCount_beta, rightCount_beta, minSamples_beta)
beta_df <- beta_df[rownames(beta_df) %in% CpGs_with_enough_coverage,]
gc()

coord_df <- as.data.frame(coord)

coord_df$end <- coord_df$start + 1
coord_df$strand <- NULL
coord_df$width <- NULL
coord_df$CpG_ID <- paste0("CpG_", 1:nrow(coord_df))
coord_df <- coord_df[coord_df$CpG_ID %in% rownames(beta_df),]

write.table(beta_df,
            file=output_path_beta,
            row.names=T,
            col.names=T,
            quote=F,
            sep="\t")

write.table(coord_df,
            file=output_path_CpGs,
            row.names=F,
            col.names=F,
            quote=F,
            sep="\t")