options(scipen=999)

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

input_path_CpGs               <- args$input_path_CpGs
output_path_ff_selected_TxIDs <- args$output_path_ff_selected_TxIDs
output_path_filtered_CpGs     <- args$output_path_filtered_CpGs
output_path_num_CpGs_plot     <- args$output_path_num_CpGs_plot
output_path_dist_filt_plot    <- args$output_path_dist_filt_plot

`%>%` <- magrittr::`%>%`

# ---- load ----
top_K_CpGs <- data.table::fread(input_path_CpGs, header = FALSE)
colnames(top_K_CpGs) <- c("chr", "start", "end", "CpGID", "TxID", "distance")

# ---- count CpGs per TxID ----
top_K_CpGs <- top_K_CpGs %>%
    dplyr::group_by(TxID) %>%
    dplyr::mutate(numCpGs = dplyr::n()) %>%
    dplyr::ungroup()

# ---- save CpGs and TxIDs ----
CpGs_filtered <- top_K_CpGs[, c("TxID", "CpGID")]

write.table(x         = CpGs_filtered,
            file      = output_path_filtered_CpGs,
            col.names = FALSE,
            row.names = FALSE,
            quote     = FALSE,
            sep       = "\t")

TxIDs <- unique(top_K_CpGs$TxID)
writeLines(text = TxIDs, sep = "\n", con = output_path_ff_selected_TxIDs)

# ---- plot: CpGs per TxID ----
pdf(output_path_num_CpGs_plot, width = 6, height = 5)

x <- unique(top_K_CpGs[, c("TxID", "numCpGs")])$numCpGs

hist(
    x,
    breaks = 30,
    col    = "steelblue",
    border = "white",
    main   = "Distribution of CpG counts per TxID",
    xlab   = "Number of CpGs",
    ylab   = "Frequency"
)
rug(x, col = "black")
abline(v = median(x, na.rm = TRUE), col = "red",       lwd = 2)
abline(v = mean(x,   na.rm = TRUE), col = "darkgreen", lwd = 2, lty = 2)
legend("topright",
       legend = c("Median", "Mean"),
       col    = c("red", "darkgreen"),
       lwd    = c(2, 2),
       lty    = c(1, 2),
       bty    = "n")

dev.off()

# ---- plot: distance from TSS ----
pdf(output_path_dist_filt_plot, width = 6, height = 5)

dist_vals <- log10(top_K_CpGs$distance)

hist(
    dist_vals,
    breaks = 45,
    col    = "darkorange",
    border = "white",
    main   = "Distribution of log10(distance from TSS)",
    xlab   = "Distance (bp)",
    ylab   = "Frequency"
)
rug(dist_vals, col = "black")
abline(v = median(dist_vals, na.rm = TRUE), col = "red",       lwd = 2)
abline(v = mean(dist_vals,   na.rm = TRUE), col = "darkgreen", lwd = 2, lty = 2)
legend("topright",
       legend = c("Median", "Mean"),
       col    = c("red", "darkgreen"),
       lwd    = c(2, 2),
       lty    = c(1, 2),
       bty    = "n")

dev.off()
