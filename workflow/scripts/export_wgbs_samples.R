suppressPackageStartupMessages({
  library(HDF5Array)
  library(SummarizedExperiment)
})

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

cpgea_path <- args$cpgea
mcrpc_path <- args$mcrpc
out_path   <- args$output_samples

cpgea <- HDF5Array::loadHDF5SummarizedExperiment(cpgea_path)
mcrpc <- HDF5Array::loadHDF5SummarizedExperiment(mcrpc_path)

s_cpgea <- colnames(cpgea)
s_mcrpc <- colnames(mcrpc)
samples <- c(s_cpgea, s_mcrpc)

dup <- samples[duplicated(samples)]
if (length(dup)) {
  stop(sprintf("Duplicate sample IDs across CPGEA/MCRPC: %s",
               paste(head(dup, 10), collapse = ", ")))
}

cat(sprintf("CPGEA samples: %d\n", length(s_cpgea)))
cat(sprintf("MCRPC samples: %d\n", length(s_mcrpc)))
cat(sprintf("Total WGBS   : %d\n", length(samples)))

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
writeLines(samples, con = out_path)
cat(sprintf("Written to %s\n", out_path))