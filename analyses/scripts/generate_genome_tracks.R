options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path_coords          <- args$input_path_coords
input_path_feat_imp        <- args$input_path_feat_imp
importance                 <- args$importance
output_path                <- args$output_path

cond <- ((length(readLines(input_path_coords)) == 0)|(length(readLines(input_path_feat_imp)) == 0))

if(cond){
    writeLines(text = "", con = output_path, sep = "")
    quit()
}

coords <- read.delim(file = input_path_coords, header = F, col.names=c("chr", "start", "end", "CpGID"))
feat_imp <- read.delim(file = input_path_feat_imp, header = T)
colnames(feat_imp)[1] <- "CpGID"

coords_with_feat_imp <- merge(coords, feat_imp, by = "CpGID")

to.return <- coords_with_feat_imp[,c("chr", "start", "end", importance)]

write.table(to.return, 
            output_path, 
            row.names = F, 
            col.names = F,
            quote = F,
            sep = "\t")
