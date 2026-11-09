# Rebuild need_gene_exp_<COHORT>.csv and gene_list_<COHORT>.txt from <COHORT>.data.Rdata with ALL expressed genes
# (no promoter-probe filter), for cohorts preprocessed with the earlier version of 02_preprocess.R.
# Usage: Rscript rebuild_gene_list.R --cohort=BRCA --work_dir=<WORK_DIR>/BRCA
a <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
setwd(a$work_dir); study <- a$cohort
load(paste0(study, ".data.Rdata"))            # data$gene.exp: genes x samples, all expressed genes after the synonym search
old <- if (file.exists(paste0("gene_list_", study, ".txt"))) length(readLines(paste0("gene_list_", study, ".txt"))) else NA
g <- unique(rownames(data$gene.exp))
write.csv(data$gene.exp[rownames(data$gene.exp) %in% g, , drop = FALSE], paste0("need_gene_exp_", study, ".csv"), row.names = TRUE)
writeLines(g, paste0("gene_list_", study, ".txt"))
cat(sprintf("%s: gene list %s -> %d genes (all expressed genes)\n", study, old, length(g)))
