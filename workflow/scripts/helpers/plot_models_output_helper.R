plot_models_output <- function(df, output_path, what){

`%>%` <- magrittr::`%>%`

numeric_df <- df %>% dplyr::select(dplyr::where(is.numeric))

df_long <- numeric_df %>%
    tidyr::pivot_longer(cols = dplyr::everything(),
                        names_to = "parameter",
                        values_to = "value")

p <- ggplot2::ggplot(df_long, ggplot2::aes(x = value)) +
    ggplot2::geom_density(fill = "#69b3a2", alpha = 0.6, color = "black") +
    ggplot2::facet_wrap(~ parameter, scales = "free", ncol = 3) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::labs(
      title = paste0("Distribution of ", what ," across all Txs"),
      x = NULL,
      y = "Density"
    )+
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

ggplot2::ggsave(
    filename = output_path,
    plot = p,
    width = 12, height = 8, dpi = 300
  )

}