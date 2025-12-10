options(repos = c(CRAN = "https://cloud.r-project.org"))
options(scipen = 999)

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")

source("workflow/scripts/helpers/plot_models_output_helpers.R")

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path_hyper         <- args$input_path_hyper
input_path_perf          <- args$input_path_perf

output_path_hyper_density        <- args$output_path_hyper_density
output_path_perf_density         <- args$output_path_perf_density
output_path_hyper_heatmap        <- args$output_path_hyper_heatmap
output_path_perf_heatmap         <- args$output_path_perf_heatmap

hypers <- read.delim(input_path_hyper, header = T)
perfs <- read.delim(input_path_perf, header = T)

plot_models_output_density(hypers, output_path_hyper_density, "Hyperparameters")
plot_models_output_heatmap(hypers, output_path_hyper_heatmap, "Hyperparameters")
plot_models_output_density(perfs, output_path_perf_density, "Performances")
plot_models_output_heatmap(perfs, output_path_perf_heatmap, "Performances")

