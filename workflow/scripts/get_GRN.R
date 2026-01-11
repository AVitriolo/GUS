options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path_composite_selected         <- args$input_path_composite_selected
input_path_pairs                      <- args$input_path_pairs
output_path_GRN                       <- args$output_path_GRN

# load files with couples and coordinates
# merge by id
cpgGenes <- read.table(input_path_pairs,head=F,sep="\t")
names(cpgGenes) <- c("id","gene")

cpgCoord <- read.table(input_path_composite_selected,head=F,sep="\t")
names(cpgCoord) <- c("chr","start","end","id")

CpGGenes <- merge(cpgGenes,cpgCoord, by="id")

write.table(CpGGenes[,c("chr","start","end","id","gene")],
                row.names=F,
                col.names=F,
                quote=F,
                sep="\t",
                file=output_path_GRN)
