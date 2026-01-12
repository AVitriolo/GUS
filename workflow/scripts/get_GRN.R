options(scipen=999)

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

input_path_gimme                     <- args$input_path_gimme
input_path_converter                 <- args$input_path_converter
input_path_slopped                   <- args$input_path_slopped
input_path_gff                       <- args$input_path_gff
input_path_expressed_TxIDs           <- args$input_path_expressed_TxIDs
output_path                          <- args$output_path

gimme <- read.table(input_path_gimme,head=F)
names(gimme) <- c("Motif","CpG","FPR")

converter <- read.table(input_path_converter,head=T,sep="\t")

gimme_with_conversion <- merge(gimme,converter,by="Motif")
gimme_with_conversion$TF <- toupper(gimme_with_conversion$Factor)

# load gene - cpg association
slopped <- read.table(input_path_slopped,head=F)
slopped <- slopped[,c(4:5)]
names(slopped) <- c("CpG","Gene")

gimme_with_conversion_with_slopped <- merge(gimme_with_conversion,slopped,by="CpG")
names(gimme_with_conversion_with_slopped) <- c("CpG","Motif","FPR","Factor","Evidence","Curated","TF","transcript_id")

gff3 <- rtracklayer::import.gff3(input_path_gff)
gff3 <- gff3[,c("gene_name","transcript_id")]
gff3$transcript_id <- gsub("\\.[0-9]*","",gff3$transcript_id)
gff3.df <- mcols(gff3)
gff3.df <- as.data.frame(gff3.df)
gff3.df.xs <- unique(gff3.df)

gimme_with_conversion_with_slopped <- gimme_with_conversion_with_slopped[,c("TF","transcript_id","FPR")]
gimme_with_conversion_with_slopped <- unique(gimme_with_conversion_with_slopped)

gimme_with_conversion_with_slopped_with_gff <- merge(gimme_with_conversion_with_slopped,gff3.df.xs,by="transcript_id")
gimme_with_conversion_with_slopped_with_gff$transcript_id <- NULL
gimme_with_conversion_with_slopped_with_gff <- (unique(gimme_with_conversion_with_slopped_with_gff))

expressed <- as.character(read.table(input_path_expressed_TxIDs,head=F)$V1)

gff3.df.xs.expr <- gff3.df.xs[which(gff3.df.xs$transcript_id %in% expressed),]

GRN <- gimme_with_conversion_with_slopped_with_gff[which(gimme_with_conversion_with_slopped_with_gff$TF %in% gff3.df.xs.expr$gene_name),]
GRN.f <- (GRN[which(abs(GRN$FPR) >= 7),])
GRN.f <- GRN.f[order(GRN.f$FPR,decreasing=T),]

GRN.ff <- GRN.f %>% 
                dplyr::group_by(TF,gene_name) %>% 
                dplyr::filter(abs(FPR) == max(abs(FPR))) %>% 
                dplyr::arrange(TF,gene_name,FPR)

GRN.fff <- GRN.ff[which(abs(GRN.ff$FPR) >= 12),]

write.table(GRN.fff,
            quote=F,
            sep="\t",
            file=output_path)
