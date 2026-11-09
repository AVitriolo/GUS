options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("limma")

options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path                               <- args$input_path
output_distance_filtered_path            <- args$output_distance_filtered_path
output_plot_1_path                       <- args$output_plot_1_path
output_plot_2_path                       <- args$output_plot_2_path
distance_cutoff                          <- as.integer(args$distance_cutoff)
min_CpG                                  <- as.integer(args$min_CpG)

CpG_Tx_distance <- read.table(input_path, header=F)
names(CpG_Tx_distance) <- c("CpGID","TxID","distance")

CpG_Tx_distance.f <- CpG_Tx_distance[which(CpG_Tx_distance$distance <= distance_cutoff),] # we could go as low as 250k

negative_distances_idxs <- which(CpG_Tx_distance.f$distance == -1)                        # we remove BedTools uncertainties

if(length(negative_distances_idxs) > 0){
    CpG_Tx_distance.ff <- CpG_Tx_distance.f[-negative_distances_idxs,]
} else {
    CpG_Tx_distance.ff <- CpG_Tx_distance.f
}
       
nCpG <- table(CpG_Tx_distance.ff$TxID)
nCpG.f <- nCpG[which(nCpG >= min_CpG)]

CpG_Tx_distance.fff <- CpG_Tx_distance.ff[which(CpG_Tx_distance.ff$TxID %in% names(nCpG.f)),]

pdf(output_plot_1_path,width=6,height=4)
limma::plotDensities(log10(nCpG.f + 1))
dev.off()

pdf(output_plot_2_path,width=6,height=4)
limma::plotDensities(log10(CpG_Tx_distance.fff$distance + 1))
dev.off()

write.table(CpG_Tx_distance.fff, 
            file=output_distance_filtered_path,
            col.names=F,
            row.names=F,
            quote=F, 
            sep="\t")

