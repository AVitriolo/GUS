options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path                               <- args$input_path
output_distance_filtered_path            <- args$output_distance_filtered_path
output_plot_1_path                       <- args$output_plot_1_path
output_plot_2_path                       <- args$output_plot_2_path

CpG_Tx_distance <- read.table(input_path, header=F)
# txids <- unique(CpG_Tx_distance$V2)

CpG_Tx_distance.f <- CpG_Tx_distance[which(CpG_Tx_distance$V3 <= 1000000),]
CpG_Tx_distance.ff <- CpG_Tx_distance.f[-which(CpG_Tx_distance.f$V3 == -1),]

nCpG <- table(CpG_Tx_distance.ff$V2)
nCpG.f <- nCpG[which(nCpG > summary(nCpG)[[2]])]

CpG_Tx_distance.fff <- CpG_Tx_distance.ff[which(CpG_Tx_distance.ff$V2 %in% names(nCpG.f)),]

pdf(output_plot_1_path,width=2,height=4)
limma::plotDensities(log10(nCpG.f))
dev.off()

pdf(output_plot_2_path,width=5,height=3)
limma::plotDensities(log10(CpG_Tx_distance.fff$V3))
dev.off()

write.table(CpG_Tx_distance.fff, 
            file=output_distance_filtered_path,
            col.names=F,
            row.names=F,
            quote=F, 
            sep="\t")

