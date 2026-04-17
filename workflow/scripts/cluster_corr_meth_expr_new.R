source("workflow/scripts/helpers/cluster_corr_meth_expr_helper_new.R")

options(scipen = 999)

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

input_path_corr  <- args$input_path_corr
output_path_corr <- args$output_path_corr
output_path_plot <- args$output_path_plot
TxID             <- args$TxID

cluster_corr(
  input_path  = input_path_corr,
  output_corr = output_path_corr,
  output_plot = output_path_plot,
  TxID        = TxID
)
