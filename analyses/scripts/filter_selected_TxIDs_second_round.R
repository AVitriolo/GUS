options(scipen=999)                                                                                            #  unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                                             #  read args

input_path_K_closest                      <- args$input_path_K_closest
distance                                  <- as.numeric(args$distance)
min_CpG                                   <- as.numeric(args$min_CpG)
output_path_ff_selected_TxIDs             <- args$output_path_ff_selected_TxIDs
output_path_filtered_CpGs                 <- args$output_path_filtered_CpGs
output_path_num_CpGs_plot                 <- args$output_path_num_CpGs_plot
output_path_dist_filt_plot                <- args$output_path_dist_filt_plot

`%>%` <- magrittr::`%>%`

#### LOAD Top K CpGs

top_K_CpGs <- data.table::fread(input_path_K_closest)
colnames(top_K_CpGs) <- c("chr", "start", "end", "CpGID", "TxID", "distance")

#### Filter top K CpGs

distance.minCpGs.filt <- top_K_CpGs %>% 
                            dplyr::filter(.data$distance <= .env$distance) %>%                                 # filter by distance .env here is needed because otherwise I compare the column with itself
                            dplyr::group_by(TxID) %>% 
                            dplyr::mutate(numCpGs = dplyr::n()) %>% 
                            dplyr::filter(numCpGs >= min_CpG)                                                  # filter by numCpGs


CpGs_filtered <- distance.minCpGs.filt[,c("TxID", "CpGID")]

write.table(x = CpGs_filtered, 
            file = output_path_filtered_CpGs,
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")

TxIDs <- unique(distance.minCpGs.filt$TxID)
writeLines(text = TxIDs, sep = "\n", con = output_path_ff_selected_TxIDs)

pdf(output_path_num_CpGs_plot, width = 6, height = 5)

x <- distance.minCpGs.filt$numCpGs

hist(
  x,
  breaks = 15,
  col = "steelblue",
  border = "white",
  main = "Distribution of CpG Counts",
  xlab = "Number of CpGs",
  ylab = "Frequency"
)

rug(x, col = "black")

abline(v = median(x, na.rm = TRUE), col = "red", lwd = 2)
abline(v = mean(x, na.rm = TRUE), col = "darkgreen", lwd = 2, lty = 2)

legend(
  "topright",
  legend = c("Median", "Mean"),
  col = c("red", "darkgreen"),
  lwd = c(2, 2),
  lty = c(1, 2),
  bty = "n"
)

dev.off()

pdf(output_path_dist_filt_plot, width = 6, height = 5)

dist_vals <- log10(distance.minCpGs.filt$distance)

hist(
  dist_vals,
  breaks = 45,
  col = "darkorange",
  border = "white",
  main = "Distribution of log10(Distance)",
  xlab = "log10(Distance)",
  ylab = "Frequency"
)

rug(dist_vals, col = "black")

abline(v = median(dist_vals, na.rm = TRUE), col = "red", lwd = 2)
abline(v = mean(dist_vals, na.rm = TRUE), col = "darkgreen", lwd = 2, lty = 2)

legend(
  "topright",
  legend = c("Median", "Mean"),
  col = c("red", "darkgreen"),
  lwd = c(2, 2),
  lty = c(1, 2),
  bty = "n"
)

dev.off()
