options(scipen=999)

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

input_dir_maps          <- args$input_dir_maps
output_path_n_clusters  <- args$output_path_n_clusters

`%>%` <- magrittr::`%>%`

# ---- find all cluster map files ----
map_files <- list.files(input_dir_maps, pattern = "_clusters_map\\.txt$", full.names = TRUE)
cat(sprintf("Found %d cluster map files\n", length(map_files)))

# ---- count distinct clusters per TxID ----
n_clusters_per_txid <- sapply(map_files, function(f) {
    mapping <- read.table(f, header = FALSE, sep = "\t")
    colnames(mapping) <- c("CpGID", "cluster_id")
    length(unique(mapping$cluster_id))
})

# ---- plot ----
pdf(output_path_n_clusters, width = 6, height = 5)

x <- n_clusters_per_txid

hist(
    x,
    breaks = 30,
    col    = "steelblue",
    border = "white",
    main   = "Distribution of clusters per TxID",
    xlab   = "Number of clusters",
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