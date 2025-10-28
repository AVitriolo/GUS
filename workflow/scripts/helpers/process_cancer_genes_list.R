options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path                               <- args$input_path
output_path 	                         <- args$output_path
assembly_code                            <- args$assembly_code

cancer_genes <- read.delim(input_path, check.names = F)
cancer_genes.f <- cancer_genes[cancer_genes[,8] > 2,]
cancer_genes.ff <- cancer_genes.f[cancer_genes.f[,7] %in% c("ONCOGENE", "ONCOGENE_AND_TSG", "TSG"),]


if(assembly_code == "hg19"){
    cancer_genes.fff <- cancer_genes.ff[,c(3)]
} else if(assembly_code == "hg38"){
    cancer_genes.fff <- cancer_genes.ff[,c(5)]
}


write.table(cancer_genes.fff, 
            file=output_path,
            col.names=F,
            row.names=F,
            quote=F, 
            sep="\t")

