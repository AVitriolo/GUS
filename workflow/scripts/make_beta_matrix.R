beta_mat <- h5_norm_prim@assays@data$beta
beta_df= as.data.frame(beta_mat)
beta_df$CpG_ID = paste0("CpG_", 1:nrow(beta_df))
write.table(beta_df, col.names = T, row.names = F, quote = F, file = "/home/Prostate_GAS/resources/20251007_CpGs_beta_NT.txt")
