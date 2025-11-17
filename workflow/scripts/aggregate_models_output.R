source("workflow/scripts/helpers/aggregate_models_output_helper.R")

options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_dir_hyperparams           <- args$input_dir_hyperparams
output_path_hyperparams         <- args$output_path_hyperparams

input_dir_performances          <- args$input_dir_performances
output_path_performances        <- args$output_path_performances

aggregate_models_output(input_dir_hyperparams, output_path_hyperparams)
aggregate_models_output(input_dir_performances, output_path_performances)

