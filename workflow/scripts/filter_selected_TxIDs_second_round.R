options(scipen=999)                                                                                            #  unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                                             #  read args

input_path_K_closest                      <- args$input_path_K_closest
distance                                  <- as.numeric(args$distance)
min_CpG                                   <- as.numeric(args$min_CpG)
output_path_ff_selected_TxIDs             <- args$output_path_ff_selected_TxIDs

`%>%` <- magrittr::`%>%`

#### LOAD Top K CpGs

top_K_CpGs <- data.table::fread(input_path_K_closest)
colnames(top_K_CpGs) <- c("chr", "start", "end", "CpGID", "TxID", "distance")

#### Filter top K CpGs

distance.minCpGs.filt <- top_K_CpGs %>% 
                            dplyr::filter(.data$distance <= .env$distance) %>%                                 # filter by distance .env here is needed because otherwise I compare the column with itself
                            dplyr::group_by(TxID) %>% 
                            dplyr::mutate(numCpGs = dplyr::n()) %>% 
                            dplyr::filter(numCpGs >= min_CpG)                                                  # filter by numCpGs

TxIDs <- unique(distance.minCpGs.filt$TxID)
writeLines(text = TxIDs, sep = "\n", con = output_path_ff_selected_TxIDs)
