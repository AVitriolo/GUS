source("workflow/scripts/helpers/plot_models_output_helper.R")

options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path_hyper         <- args$input_path_hyper
input_path_perf          <- args$input_path_perf

output_path_hyper        <- args$output_path_hyper
output_path_perf        <- args$output_path_perf

hypers <- read.delim(input_path_hyper, header = T)
perfs <- read.delim(input_path_perf, header = T)

plot_models_output(hypers, output_path_hyper, "Hyperparameters")
plot_models_output(perfs, output_path_perf, "Performances")

