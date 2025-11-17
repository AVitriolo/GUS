options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_dir          <- args$input_dir
output_path        <- args$output_path

hyper_output <- list.files(input_dir, pattern = "ENST", full.names = T)

TxIDs <- unlist(lapply(X = hyper_output, 
                       FUN = function(x){
                        parts <- strsplit(x, split = "_")[[1]] 
                        tail(parts, n = 1)})
                        )

hyperparams_list <- list()

for(idx in seq(hyper_output)){
    TxID <- TxIDs[idx]
    file <- hyper_output[idx]
    tab <- read.delim(file, header = F)
    rownames(tab) <- tab$V1
    tab$V1 <- NULL
    colnames(tab) <- TxID
    hyperparams_list[[idx]] <- as.data.frame(t(tab))
}

hyperparams_df <- do.call(rbind, hyperparams_list)
hyperparams_df$TxID <- rownames(hyperparams_df)
rownames(hyperparams_df) <- NULL

hyperparams_df <- hyperparams_df[, c(ncol(hyperparams_df), 1:(ncol(hyperparams_df) - 1))]

write.table(hyperparams_df, 
            output_path, 
            row.names = F, 
            col.names = T,
            quote = F,
            sep = "\t")
