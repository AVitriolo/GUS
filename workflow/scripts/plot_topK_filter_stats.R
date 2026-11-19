# QC plots of the candidate-CpG sets, from the aggregates of split_topK_by_TxID.sh
# (same two figures as filter_selected_TxIDs_second_round_NTM.R, without loading the CpG table).
# Usage: Rscript plot_topK_filter_stats.R <cpgs_per_txid.tsv> <distance_hist.tsv> <num_CpGs.pdf> <distance.pdf>
args <- commandArgs(trailingOnly = TRUE)
counts <- read.table(args[1], sep = "\t", col.names = c("TxID", "numCpGs"))
dhist  <- read.table(args[2], sep = "\t", col.names = c("log10_bin", "n"))
dhist  <- dhist[order(dhist$log10_bin), ]

pdf(args[3], width = 6, height = 5)
x <- counts$numCpGs
hist(x, breaks = 30, col = "steelblue", border = "white",
     main = "Distribution of CpG counts per TxID", xlab = "Number of CpGs", ylab = "Frequency")
rug(x, col = "black")
abline(v = median(x), col = "red", lwd = 2)
abline(v = mean(x), col = "darkgreen", lwd = 2, lty = 2)
legend("topright", legend = c("Median", "Mean"), col = c("red", "darkgreen"), lwd = 2, lty = c(1, 2), bty = "n")
dev.off()

pdf(args[4], width = 6, height = 5)
mids <- dhist$log10_bin + 0.05
wmean <- sum(mids * dhist$n) / sum(dhist$n)
cum <- cumsum(dhist$n); wmed <- mids[which(cum >= sum(dhist$n) / 2)[1]]
barplot(dhist$n, names.arg = sprintf("%.1f", dhist$log10_bin), col = "darkorange", border = "white", space = 0,
        main = "Distribution of log10(distance from TSS)", xlab = "log10 distance (bp), 0.1 bins", ylab = "CpG-TSS pairs",
        las = 2, cex.names = 0.6)
abline(v = (wmed - min(dhist$log10_bin)) * 10, col = "red", lwd = 2)
abline(v = (wmean - min(dhist$log10_bin)) * 10, col = "darkgreen", lwd = 2, lty = 2)
legend("topleft", legend = c("Median", "Mean"), col = c("red", "darkgreen"), lwd = 2, lty = c(1, 2), bty = "n")
dev.off()
cat(sprintf("%d TxIDs; CpGs per TxID: median %.0f, mean %.0f\n", nrow(counts), median(x), mean(x)))
