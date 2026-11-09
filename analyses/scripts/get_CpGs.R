options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path             <- args$input_path
output_path            <- args$output_path

rse <- HDF5Array::loadHDF5SummarizedExperiment(input_path)

coord  <- SummarizedExperiment::rowRanges(rse)
coord_df <- as.data.frame(coord)

coord_df$end <- coord_df$start + 1
coord_df$strand <- NULL
coord_df$width <- NULL
coord_df$CpG_ID <- paste0("CpG_", 1:nrow(coord_df))
 
write.table(coord_df,
            file=output_path,
            row.names=F,
            col.names=F,
            quote=F,
            sep="\t")