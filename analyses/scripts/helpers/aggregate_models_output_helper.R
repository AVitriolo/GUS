aggregate_models_output <- function(input_dir, output_path){

output <- list.files(input_dir, pattern = "ENST", full.names = T)

TxIDs <- unlist(lapply(X = output, 
                       FUN = function(x){
                        parts <- strsplit(x, split = "_")[[1]] 
                        tail(parts, n = 1)})
                        )

list <- list()

for(idx in seq(output)){
    TxID <- TxIDs[idx]
    file <- output[idx]
    tab <- read.delim(file, header = F)
    rownames(tab) <- tab$V1
    tab$V1 <- NULL
    colnames(tab) <- TxID
    tab <- as.data.frame(t(tab))
    tab$TxID <- rownames(tab)
    rownames(tab) <- NULL
    list[[idx]] <- tab
}

df <- do.call(rbind, list)

df <- df[, c(ncol(df), 1:(ncol(df) - 1))]

write.table(df, 
            output_path, 
            row.names = F, 
            col.names = T,
            quote = F,
            sep = "\t")

}
