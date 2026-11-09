# List the TCGA cohorts on the UCSC Xena tcgaHub that have >= MIN samples with BOTH HumanMethylation450 and HiSeqV2.
# Writes <out>: one cohort code per line (e.g. BRCA), and <out>.tsv with the counts for every cohort.
# Usage: Rscript 00_select_cohorts.R --min_samples=100 --out=cohorts_selected.txt
suppressPackageStartupMessages({library(UCSCXenaTools); library(dplyr)})
a <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
MIN <- as.integer(ifelse(is.null(a$min_samples), 100, a$min_samples)); out <- ifelse(is.null(a$out), "cohorts_selected.txt", a$out)
host <- "https://tcga.xenahubs.net"
d <- XenaData %>% filter(XenaHostNames == "tcgaHub", grepl("^TCGA\\.[A-Z]+\\.sampleMap/(HumanMethylation450|HiSeqV2)$", XenaDatasets))
d$code <- sub("^TCGA\\.([A-Z]+)\\.sampleMap/.*", "\\1", d$XenaDatasets)
codes <- sort(unique(d$code[grepl("HumanMethylation450$", d$XenaDatasets)]))
codes <- codes[codes %in% d$code[grepl("HiSeqV2$", d$XenaDatasets)]]
cat(sprintf("[%s] %d cohorts have both datasets on Xena; counting samples...\n", format(Sys.time(), "%H:%M:%S"), length(codes))); flush.console()
res <- do.call(rbind, lapply(seq_along(codes), function(i) { cc <- codes[i]
  s_m <- tryCatch(fetch_dataset_samples(host, paste0("TCGA.", cc, ".sampleMap/HumanMethylation450")), error = function(e) character(0))
  s_e <- tryCatch(fetch_dataset_samples(host, paste0("TCGA.", cc, ".sampleMap/HiSeqV2")), error = function(e) character(0))
  nb <- length(intersect(s_m, s_e))
  cat(sprintf("[%s] %2d/%d %-6s methylation %4d | expression %4d | both %4d%s\n", format(Sys.time(), "%H:%M:%S"), i, length(codes), cc, length(s_m), length(s_e), nb,
              ifelse(length(s_m) == 0 || length(s_e) == 0, "  (sample list could not be fetched)", ""))); flush.console()
  data.frame(cohort = cc, n_methylation = length(s_m), n_expression = length(s_e), n_both = nb)
}))
res <- res[order(-res$n_both), ]
write.table(res, paste0(out, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
# pan-cancer / combined cohorts (e.g. PANCAN, COADREAD, GBMLGG, LUNG) overlap the single cohorts: kept in the .tsv, not selected
combined <- c("PANCAN", "PANCAN12", "COADREAD", "GBMLGG", "LUNG", "FPPP")
sel <- res$cohort[res$n_both >= MIN & !res$cohort %in% combined]
writeLines(sel, out)
print(res, row.names = FALSE)
cat(sprintf("\n%d cohorts with >= %d samples having both data types: %s\n", length(sel), MIN, paste(sel, collapse = " ")))
