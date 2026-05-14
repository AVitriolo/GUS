source("workflow/scripts/helpers/find_neighbors.R")
source("workflow/scripts/helpers/connect_KNN_igraph.R")
source("workflow/scripts/helpers/cluster_CpGs.R")

set.seed(123)

options(scipen=999)                                                                 # unable scientific notation           

library(FNN)
library(ggplot2)

#args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                 # read args

input_path_xgb_input           <- "./xgb_hg38_T_0.5_0.5_90_10_v38_pc_20000_50_25_10_1000_ENST00000003084"
input_path_dataset_corr        <- "./corr_hg38_T_0.5_0.5_90_10_v38_pc_20000_50_25_10_1000_ENST00000003084"
TxID			               <- "ENST00000003084"
C                              <- 2
B                              <- 10

corr <- read.delim(input_path_dataset_corr, check.names = F, header = T, sep = "\t") # corr matrix
corr <- abs(corr)

xgb <- read.delim(input_path_xgb_input, check.names = F, header = T, sep = "\t") # xgb_input
xgb <- xgb[,-ncol(xgb)] # remove TxID

percentiles_list <- seq(from = 10, to = 1000, length.out = 100)       # define ranges
flex_points_list <- seq(from = 1, to = 3, by = 1)                     # define ranges
params_grid <- expand.grid(percentiles = percentiles_list, flex_points = flex_points_list)

CpGIDs <- colnames(corr)
chrs <- unlist(lapply(X = CpGIDs, FUN = function(x){strsplit(x = x, split = "_")[[1]][1]}))
starts <- unlist(lapply(X = CpGIDs, FUN = function(x){strsplit(x = x, split = "_")[[1]][2]}))
coords <- data.frame(chr = chrs, start = starts, end = starts, CpGID = CpGIDs)

res <- cluster_CpGs(params_grid, corr, coords, C)
n_clusters_observed <- res$n_clusters

n_clusters_expected_vec <- parallel::mclapply(
    X = 1:B,
    mc.cores = C,
    FUN = function(b){

        print(b)

        colnames(xgb) <- sample(colnames(xgb))

        random_corr <- abs(cor(xgb, method = "spearman"))

        num_part <- as.numeric(sub(".*_(\\d+)$", "\\1", colnames(random_corr)))

        random_corr <- random_corr[order(num_part), order(num_part)]

        stopifnot(rownames(random_corr) == rownames(corr))

        random_CpGIDs <- colnames(random_corr)

        random_coords <- data.frame(
            chr = sapply(strsplit(random_CpGIDs, "_"), `[`, 1),
            start = sapply(strsplit(random_CpGIDs, "_"), `[`, 2),
            end = sapply(strsplit(random_CpGIDs, "_"), `[`, 2),
            CpGID = random_CpGIDs
        )

        random_res <- cluster_CpGs(params_grid, random_corr, random_coords, C)

        return(random_res$n_clusters)
    }
)

n_clusters_expected_vec <- unlist(n_clusters_expected_vec)

pseudo_pval <- sum(n_clusters_expected_vec <= n_clusters_observed) / B

histo <- ggplot(data = data.frame(n_clusters_expected_vec), aes(x = n_clusters_expected_vec)) +
geom_density(fill = "skyblue", alpha = 0.5) +
geom_vline(xintercept = n_clusters_observed, color = "red", linetype = "dashed", linewidth = 1) +
labs(title = "Expected vs Observed Number of Clusters", x = "Number of Clusters", y = "Density") +
annotate("text", x = Inf, y = Inf, label = paste0("pseudo p-value: ", pseudo_pval), hjust = 1.1, vjust = 1.5, size = 5) +
theme_bw()

ggsave(paste0(TxID, "_histo.png"), histo, width = 7, height = 7)
