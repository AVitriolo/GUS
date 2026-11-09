# Preprocess one cohort exactly as in the TCGA PRAD analysis:
#   data.processing() of read.data.preprocessing_function_infinite.r (beta -> M, beta 0/1 -> -10/10, REMP grooMethy
#   imputation of probes with <= 20% NA, drop probes still NA, keep genes with mean expression > 1, matched samples,
#   hg19 450K annotation, NCBI synonym search), then the exports of code/0_download_data.R / code/2_export_csv.R
#   (row names kept, probe-level ranges, expression not transposed) and the gene list.
#   Difference from the original: the gene list is ALL expressed genes (no requirement of a probe in the promoter).
#   The methylation matrix is written twice: M-values (<STUDY>_data_meth.csv, as before) and beta values
#   (<STUDY>_data_meth_beta.csv = 2^M/(1+2^M) of the same imputed matrix; imputation itself is done on M, unchanged).
# Only the hard-coded paths were turned into arguments; the processing itself is unchanged.
# Usage: Rscript 02_preprocess.R --cohort=BRCA --work_dir=<WORK_DIR>/BRCA --resources=<PIPE_DIR>/resources --scripts=<PIPE_DIR>/scripts
suppressPackageStartupMessages({library(data.table); library(REMP); library(GenomicRanges); library(parallel)})
a <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
study <- a$cohort; RES <- a$resources; ncores <- as.integer(ifelse(is.null(a$cores), 8, a$cores))
setwd(a$work_dir)
GENE_INFO_PATH <- file.path(RES, "Homo_sapiens.gene_info")
PROMOTER_PATH <- file.path(RES, "hg19_promoter.txt")
source(file.path(a$scripts, "synonym.search.r"))

data.processing <- function(study) {
  f450 <- paste0("TCGA.", study, ".sampleMap/HumanMethylation450_", study)
  meth.beta <- fread(f450, header = TRUE, check.names = FALSE)
  sj.names <- read.table(f450, header = FALSE, row.names = 1, check.names = FALSE, nrows = 1, stringsAsFactors = FALSE)[1, ]
  meth.beta <- data.frame(meth.beta); rownames(meth.beta) <- meth.beta[, 1]; meth.beta <- meth.beta[, -1]; colnames(meth.beta) <- sj.names
  meth.m <- log2(meth.beta / (1 - meth.beta))
  meth.m[meth.beta == 1] <- 10        # remove infinite cases
  meth.m[meth.beta == 0] <- -10
  n <- ncol(meth.m)
  percent.na <- apply(meth.m, 1, function(x) sum(is.na(x)) / n)
  impute.ind <- percent.na != 0 & percent.na <= 0.2
  if (sum(impute.ind) > 30000) {
    temp <- grooMethy(meth.m[impute.ind, ]); imputed <- getM(temp)
    meth.m.imputed <- meth.m; meth.m.imputed[impute.ind, ] <- imputed
  } else {
    impute.ind2 <- percent.na <= 0.2
    temp <- grooMethy(meth.m[impute.ind2, ]); imputed <- getM(temp)
    meth.m.imputed <- meth.m; meth.m.imputed[impute.ind2, ] <- imputed
  }
  any.na <- apply(meth.m.imputed, 1, function(x) any(is.na(x)))
  cat(sprintf("probes with NA after imputation (dropped): %.1f%%\n", 100 * sum(any.na) / nrow(meth.m.imputed)))
  meth.m.imputed.no.na <- meth.m.imputed[!any.na, ]

  gene.exp0 <- read.table(paste0("TCGA.", study, ".sampleMap/HiSeqV2_", study), header = TRUE, row.names = 1, check.names = FALSE)
  gene.exp.mean <- rowMeans(gene.exp0)
  high.gene.exp <- gene.exp0[gene.exp.mean > 1, ]
  meth0 <- meth.m.imputed.no.na[, colnames(meth.m.imputed.no.na) %in% colnames(high.gene.exp)]
  gene.exp <- high.gene.exp[, colnames(meth0)]
  cat(sprintf("matched samples: %d | expressed genes: %d\n", ncol(gene.exp), nrow(gene.exp)))
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  probe <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19, what = "Locations")
  probe.range <- GRanges(seqnames = probe$chr, ranges = IRanges(probe$pos, width = 1, names = rownames(probe)), strand = probe$strand, name = rownames(probe))
  id <- match(rownames(meth0), names(probe.range)); sort.non.na.id <- sort(id[!is.na(id)])
  meth.range <- probe.range[sort.non.na.id]
  meth <- meth0[names(meth.range), ]
  data <- list(meth = meth, gene.exp = gene.exp, meth.range = meth.range)

  promoter <- read.delim(PROMOTER_PATH)
  promoter.range <- GRanges(seqnames = promoter$chrID, ranges = IRanges(start = promoter$start, end = promoter$end), strand = promoter$strand, gene.name = promoter$gene.name)
  syn.name <- synonym.search(promoter.range$gene.name, rownames(data$gene.exp))
  data$gene.exp <- as.matrix(data$gene.exp); rownames(data$gene.exp) <- syn.name
  save(data, file = paste0(study, ".data.Rdata"))

  probe.range <- data$meth.range
  ids <- unique(subjectHits(findOverlaps(probe.range, promoter.range))); sids <- sort(ids)
  probe.per.promoter <- mclapply(seq_along(sids), function(curi) probe.range[queryHits(findOverlaps(probe.range, promoter.range[sids[curi]]))], mc.cores = ncores)
  names(probe.per.promoter) <- promoter$gene.name[sids]
  save(probe.per.promoter, file = paste0("probe.per.promoter_", study, ".Rdata"))
  list(data = data, probe.per.promoter = probe.per.promoter, promoter.range = promoter.range)
}

r <- data.processing(study)
data <- r$data; probe.per.promoter <- r$probe.per.promoter; promoter.range <- r$promoter.range
# Gene list: ALL expressed genes (mean > 1, after the synonym search). The original pipeline additionally required >= 1
# probe inside the gene's promoter (hg19_promoter.txt); that filter is NOT applied here. Genes without a TSS anchor in
# the gene-association table are skipped later, when the per-gene tables are built.
gene.list2 <- unique(rownames(data$gene.exp))
write.csv(data$meth, paste0(study, "_data_meth.csv"), row.names = TRUE)                       # probes x samples, M-values
write.csv(2^data$meth / (1 + 2^data$meth), paste0(study, "_data_meth_beta.csv"), row.names = TRUE)   # same probes x samples, beta = 2^M/(1+2^M)
meth.range.df <- as.data.frame(data$meth.range); meth.range.df$probe.id <- names(data$meth.range)
write.csv(meth.range.df, paste0(study, "_data_meth_range.csv"), row.names = TRUE)             # probe-level ranges (hg19)
need_gene_exp_all <- data$gene.exp[rownames(data$gene.exp) %in% gene.list2, , drop = FALSE]
write.csv(need_gene_exp_all, paste0("need_gene_exp_", study, ".csv"), row.names = TRUE)        # genes x samples
writeLines(unique(rownames(need_gene_exp_all)), paste0("gene_list_", study, ".txt"))
cat(sprintf("Wrote %s_data_meth.csv (%d probes x %d samples), %s_data_meth_range.csv, need_gene_exp_%s.csv (%d genes), gene_list_%s.txt\n",
            study, nrow(data$meth), ncol(data$meth), study, study, nrow(need_gene_exp_all), study))
