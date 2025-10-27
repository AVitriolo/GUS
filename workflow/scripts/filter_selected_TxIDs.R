options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path_selected_TxIDs                               <- args$input_path_selected_TxIDs
input_path_filtered_TxIDs_num_CpGs                      <- args$input_path_filtered_TxIDs_num_CpGs
input_path_kallisto                                     <- args$input_path_kallisto
output_path                                             <- args$output_path


selected_TxIDs <- readLines(input_path_selected_TxIDs)

filtered_TxIDs_num_CpGs <- read.table(input_path_filtered_TxIDs_num_CpGs, header=F)
names(filtered_TxIDs_num_CpGs) <- c("CpG_ID","TxID","distance")
filtered_TxIDs_num_CpGs <- filtered_TxIDs_num_CpGs$TxID

counts <- read.table(input_path_kallisto,head=T)
kallisto_TxIDs <- rownames(counts)

if (length(selected_TxIDs) == 0){

    TxIDs_to_process <- intersect(filtered_TxIDs_num_CpGs, kallisto_TxIDs)

} else {

    TxIDs_to_process <- selected_TxIDs[selected_TxIDs %in% intersect(filtered_TxIDs_num_CpGs, kallisto_TxIDs)]
    TxIDs_to_discard <- selected_TxIDs[!(selected_TxIDs %in% intersect(filtered_TxIDs_num_CpGs, kallisto_TxIDs))]

    for (TxID in TxIDs_to_discard){
        msg <- paste0("Selected TxID ", TxID, " either is missing from RNA counts or does not have minCpG CpGs, skipping its execution")
        message(msg)
    }

}

writeLines(text = TxIDs_to_process, con = output_path, sep = "\n")