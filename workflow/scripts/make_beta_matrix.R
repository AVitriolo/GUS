options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path             <- args$input_path
output_path            <- args$output_path
sample_type            <- args$sample_type

h5_norm_prim <- methrix::load_HDF5_methrix(input_path)

beta_mat <- h5_norm_prim@assays@data$beta
beta_df <- as.data.frame(beta_mat)

if(sample_type == "T"){

beta_df <- beta_df[, grep("^T", colnames(beta_df))]
beta_df$CpG_ID = paste0("CpG_", 1:nrow(beta_df))

} else if (sample_type == "N"){

beta_df <- beta_df[, grep("^N", colnames(beta_df))]
beta_df$CpG_ID = paste0("CpG_", 1:nrow(beta_df))

}

write.table(beta_df,
            file=output_path,
            row.names=F,
            col.names=T,
            quote=F,
            sep="\t")
