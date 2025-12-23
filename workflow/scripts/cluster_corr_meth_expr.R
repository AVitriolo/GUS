source("workflow/scripts/helpers/cluster_corr_meth_expr_helper.R")

options(scipen = 999)

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

input_path_corr_pos  <- args$input_path_corr_pos
input_path_corr_neg  <- args$input_path_corr_neg
output_path_corr_pos <- args$output_path_corr_pos
output_path_corr_neg <- args$output_path_corr_neg
output_path_plot_pos <- args$output_path_plot_pos
output_path_plot_neg <- args$output_path_plot_neg
TxID                 <- args$TxID


cluster_corr(
  input_path  = input_path_corr_pos,
  output_corr = output_path_corr_pos,
  output_plot = output_path_plot_pos,
  TxID        = TxID,
  sign        = "pos"
)


cluster_corr(
  input_path  = input_path_corr_neg,
  output_corr = output_path_corr_neg,
  output_plot = output_path_plot_neg,
  TxID        = TxID,
  sign        = "neg"
)