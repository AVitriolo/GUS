 options(scipen = 999)
 library(methrix)
 h5_norm_prim <- load_HDF5_methrix("./cpgea_wgbs_with_coverage_hg38")
 coord = rowData(h5_norm_prim)
 coord_df = as.data.frame(coord)
 coord_df$end = coord_df$start + 1
 coord_df$strand = NULL
 coord_df$CpG_ID = paste0("CpG_", 1:nrow(coord_df))
 write.table(coord_df, row.names = F, col.names = F, sep = '\t', file = "/home/Prostate_GAS/20251007_CpGs_coord_NT.bed", quote = F)
