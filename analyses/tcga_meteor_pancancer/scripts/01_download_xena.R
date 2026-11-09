# Download TCGA HumanMethylation450 and HiSeqV2 (gene-level RSEM, log2(x+1)) for one cohort from the UCSC Xena tcgaHub.
# Same queries as the TCGA PRAD analysis (code/0_download_data.R), with the cohort as an argument.
# Usage: Rscript 01_download_xena.R --cohort=BRCA --out_dir=<WORK_DIR>/BRCA
suppressPackageStartupMessages({library(UCSCXenaTools); library(dplyr); library(R.utils)})
a <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
study <- a$cohort; out <- a$out_dir
dir.create(out, recursive = TRUE, showWarnings = FALSE); setwd(out)
options(timeout = 3600)

cohort <- unique(XenaData$XenaCohorts[XenaData$XenaHostNames == "tcgaHub" & grepl(paste0("\\(", study, "\\)$"), XenaData$XenaCohorts)])
if (length(cohort) != 1) stop(sprintf("could not find a unique tcgaHub cohort for '%s' (found: %s)", study, paste(cohort, collapse = ", ")))
cat("Xena cohort:", cohort, "\n")

for (ds in c("HumanMethylation450", "HiSeqV2$")) {
  q <- XenaGenerate(subset = XenaHostNames == "tcgaHub" & XenaCohorts == cohort) %>% XenaFilter(filterDatasets = ds)
  if (nrow(q@datasets) == 0) stop(sprintf("dataset %s not available for %s", ds, study))
  q %>% XenaQuery() %>% XenaDownload(destdir = ".")
}
smap <- paste0("TCGA.", study, ".sampleMap")
# the preprocessing function expects <file>_<STUDY> (as in TCGA.PRAD.sampleMap/HumanMethylation450_PRAD)
for (f in c("HumanMethylation450", "HiSeqV2")) {
  gunzip(file.path(smap, paste0(f, ".gz")), destname = file.path(smap, paste0(f, "_", study)), remove = FALSE, overwrite = TRUE)
}
cat("Downloaded and unpacked:", list.files(smap), "\n")
