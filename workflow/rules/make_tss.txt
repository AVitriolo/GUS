system("wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gff3.gz")
# Read in gencode.v38.annotation.gff3.gz as a GRanges 
system.time({gencode_gr = rtracklayer::import.gff3("gencode.v49.annotation.gff3.gz")})
mcols(gencode_gr) = mcols(gencode_gr)[c("ID", "gene_id", "gene_name", "transcript_support_level", "tag", "type", "transcript_type")]

# Convert tag column to a character. Took 4 minutes. 
system.time({gencode_gr$tag = sapply(gencode_gr$tag, function(x) paste(x, collapse = "; "))})

# Remove transcript version but keep PAR_Y designation and add PAR_Y to gene name
gencode_gr$ID = gsub("\\.[0-9]*", "", gencode_gr$ID)
names(gencode_gr) = gencode_gr$ID

# Create GRanges objects for protein-coding and transcripts (excluding mitochondrial transcripts)
pc_transcripts_gr = gencode_gr[gencode_gr$type == "transcript" & 
                                 gencode_gr$transcript_type == "protein_coding" & seqnames(gencode_gr) != "chrM"]

write.table(cbind(as.data.frame(pc_transcripts_tss_gr@seqnames),as.data.frame(pc_transcripts_tss_gr@ranges))[,c(1:3,5)],row.names=F,col.names=F,quote=F,file="/home/ReferenceData/20251007_hg38_gencode_v49_tss_pc.bed",sep="\t")
