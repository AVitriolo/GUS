options(scipen=999)

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

input_path_corr                           <- args$input_path_corr
input_path_sel_feat                       <- args$input_path_sel_feat
input_path_coord                          <- args$input_path_coord
output_path_cor_neg                       <- args$output_path_cor_neg
output_path_cor_pos                       <- args$output_path_cor_pos

corr = read.table(input_path_corr, header=TRUE, sep="\t")
sel_feats = read.table(input_path_sel_feat, header=FALSE, sep="\t")
coord = data.table::fread(input_path_coord, header=FALSE, sep="\t")

colnames(coord) = c("chr", "start", "end", "CpG_ID")
colnames(sel_feats) = c("CpG_ID")
colnames(corr) = c("CpG_ID", "correlation")

corr.f = corr[corr$CpG_ID %in% sel_feats$CpG_ID, ]

idx <- match(corr.f$CpG_ID, coord$CpG_ID)

corr.f$chr   <- coord$chr[idx]
corr.f$start <- coord$start[idx]
corr.f$end   <- coord$end[idx]

corr.f.pos <- corr.f[corr.f$correlation >= 0, ]
corr.f.neg <- corr.f[corr.f$correlation < 0, ]

bed.neg <- data.frame(
  chr   = corr.f.neg$chr,
  start = corr.f.neg$start,
  end   = corr.f.neg$end,
  name  = corr.f.neg$CpG_ID,
  score = corr.f.neg$correlation
)
write.table(bed.neg, file = output_path_cor_neg, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

bed.pos <- data.frame(
  chr   = corr.f.pos$chr,
  start = corr.f.pos$start,
  end   = corr.f.pos$end,
  name  = corr.f.pos$CpG_ID,
  score = corr.f.pos$correlation
)

write.table(bed.pos, file = output_path_cor_pos, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
