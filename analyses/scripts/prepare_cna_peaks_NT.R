# =============================================================================
# prepare_cna_peaks_NT.R
# CPGEA GISTIC2 focal lesions (Li et al., Nature 2020, Suppl. Table 6,
# sheet "CPGEA focal lesions") -> long table of per-sample CNA calls.
#
#   - intervals  : "Peak Limits" column (1-based, closed)
#   - call       : thresholded rows (0 = none, 1 = low-level, 2 = high-level)
#   - log2_ratio : matching "<peak> - CN values" rows (bulk log2 copy ratio)
#
# Only (sample, peak) pairs with a non-zero thresholded call are written:
# everything else is treated as copy-neutral downstream.
# =============================================================================
suppressPackageStartupMessages({ library(readxl); library(GenomicRanges) })

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

input_xlsx  <- args$input_xlsx
sheet       <- args$sheet
output_path <- args$output_path

d <- as.data.frame(suppressMessages(read_excel(input_xlsx, sheet = sheet)))
d[[1]] <- trimws(d[[1]])

sample_cols <- grep("_WGS$", colnames(d), value = TRUE)
is_cn       <- grepl(" - CN values$", d[[1]])
calls       <- d[!is_cn, ]
cn          <- d[is_cn, ]
rownames(cn) <- sub(" - CN values$", "", cn[[1]])

if (!setequal(calls[[1]], rownames(cn)))
    stop("Thresholded and CN-value rows do not describe the same peaks.")
cn <- cn[calls[[1]], ]

# "chr8:125572001-125639999(probes 158703:158709)" -> chr / start / end
limits <- sub("\\(.*$", "", calls[["Peak Limits"]])
peaks  <- GRanges(limits)
peaks$peak <- calls[[1]]
peaks$type <- ifelse(grepl("^Amplification", calls[[1]]), "amp", "del")

# a CpG must map to a single peak, otherwise its copy number is ambiguous
if (any(countOverlaps(peaks, peaks) > 1))
    stop("Overlapping Peak Limits: the per-CpG copy number would be ambiguous.")

long <- do.call(rbind, lapply(sample_cols, function(s) {
    call <- as.numeric(calls[[s]])
    keep <- call != 0
    if (!any(keep)) return(NULL)
    data.frame(
        sample     = sub("_WGS$", "", s),
        seqnames   = as.character(seqnames(peaks))[keep],
        start      = start(peaks)[keep],
        end        = end(peaks)[keep],
        peak       = peaks$peak[keep],
        type       = peaks$type[keep],
        call       = call[keep],
        log2_ratio = as.numeric(cn[[s]])[keep]
    )
}))

write.table(long, output_path, sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("%d peaks, %d samples with >=1 call, %d (sample, peak) calls -> %s\n",
            length(peaks), length(unique(long$sample)), nrow(long), output_path))
