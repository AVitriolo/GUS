options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path                               <- args$input_path
output_path                              <- args$output_path


mappability <- data.table::fread(input_path, skip = 1, header = F, data.table = F)
high_mappability <- mappability[mappability$V4 == 1, ]

write.table(high_mappability,
            file = output_path,
            row.names = F,
            col.names = F,
            quote = F,
            sep = "\t"

)