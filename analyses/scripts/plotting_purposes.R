options(scipen = 999)

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path                       <- args$input_path

output_path_nTxIDxCpG            <- args$output_path_nTxIDxCpG
output_path_nCpGxTxID            <- args$output_path_nCpGxTxID
output_path_distance             <- args$output_path_distance

surviving <- data.table::fread(input_path, header = T, sep = "\t")
surviving_unique <- dplyr::distinct(surviving)

tx_per_cpg <- surviving_unique |>
  dplyr::group_by(CpGID) |>
  dplyr::summarise(n_TxID = dplyr::n(), .groups = "drop")

cpg_per_tx <- surviving_unique |>
  dplyr::group_by(TxID) |>
  dplyr::summarise(n_CpG = dplyr::n(), .groups = "drop")

plot_tx_per_cpg <- ggplot2::ggplot(tx_per_cpg, ggplot2::aes(x = n_TxID)) +
  ggplot2::geom_histogram(binwidth = 1, boundary = 0, color = "black") +
  ggplot2::labs(
    title = "Number of TxID per CpG",
    x = "Number of TxID per CpG",
    y = "Count"
  ) +
  ggplot2::theme_minimal()

plot_cpg_per_tx <- ggplot2::ggplot(cpg_per_tx, ggplot2::aes(x = n_CpG)) +
  ggplot2::geom_histogram(binwidth = 1, boundary = 0, color = "black") +
  ggplot2::labs(
    title = "Number of CpG per TxID",
    x = "Number of CpG per TxID",
    y = "Count"
  ) +
  ggplot2::theme_minimal()

plot_distance <- ggplot2::ggplot(
  surviving_unique,
  ggplot2::aes(x = distance)
) +
  ggplot2::geom_histogram(bins = 50, color = "black") +
  ggplot2::labs(
    title = "Distribution of distance",
    x = "Distance",
    y = "Count"
  ) +
  ggplot2::theme_minimal()


ggplot2::ggsave(
  filename = output_path_nTxIDxCpG,
  plot = plot_tx_per_cpg,
  width = 8,
  height = 6,
  dpi = 300
)

ggplot2::ggsave(
  filename = output_path_nCpGxTxID,
  plot = plot_cpg_per_tx,
  width = 8,
  height = 6,
  dpi = 300
)

ggplot2::ggsave(
  filename = output_path_distance,
  plot = plot_distance,
  width = 8,
  height = 6,
  dpi = 300
)
