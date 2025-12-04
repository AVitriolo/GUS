source("workflow/scripts/helpers/load_h5_rse.R")

options(scipen=999)                                                               #  unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                #  read args

input_dir                    <- args$input_dir                                    #  resources/cpgea_wgbs_with_coverage_hg38/
input_path_mappability       <- args$input_path_mappability                       #  resources/reference_data/mappability/hg38_k100.bismap.filtered.bedgraph
output_path                  <- args$output_path

#### LOAD CpGs

h5_list <- load_h5_rse(input_dir)
rse <- h5_list$rse
CpGs <- h5_list$CpGs

#### LOAD mappability

mappability <- data.table::fread(input_path_mappability, skip = 1, header = F, data.table = F)
colnames(mappability) <- c("chr", "start", "end", "mappability")
mappability.gr <- GenomicRanges::GRanges(mappability)

#### filter CpGs by mappability

hits <- GenomicRanges::findOverlaps(query =CpGs, subject=mappability.gr)
CpGs.f <- CpGs[S4Vectors::queryHits(hits)]

CpGs.df <- as.data.frame(CpGs.f)[,c(1,2,3,6)]
CpGs.df$end <- CpGs.df$end + 1

write.table(x = CpGs.df, file = output_path, col.names = F, row.names = F, quote = F, sep = "\t")
