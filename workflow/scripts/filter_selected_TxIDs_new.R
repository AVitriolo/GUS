options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path_selected_TxIDs                               <- args$input_path_selected_TxIDs
input_path_kallisto                                     <- args$input_path_kallisto
output_path                                             <- args$output_path


selected_TxIDs <- readLines(input_path_selected_TxIDs)
counts_TxIDs <- readLines(input_path_kallisto)

if (length(selected_TxIDs) == 0){

    TxIDs_to_process <- counts_TxIDs

} else {

    TxIDs_to_process <- intersect(selected_TxIDs, counts_TxIDs)
    TxIDs_to_discard <- setdiff(selected_TxIDs, counts_TxIDs)

    for (TxID in TxIDs_to_discard){
        msg <- paste0("Selected TxID ", TxID, " either is missing from RNA counts, skipping its execution")
        message(msg)
    }

}

writeLines(text = TxIDs_to_process, con = output_path, sep = "\n")