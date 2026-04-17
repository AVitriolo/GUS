cluster_corr <- function(input_path, output_corr, output_plot, TxID) {

  plot_df  <- NULL
  write_df <- NULL

  if (!file.exists(input_path) || file.size(input_path) == 0) {

    write_df <- data.frame(
      chr         = character(),
      start       = numeric(),
      end         = numeric(),
      CpG_ID      = character(),
      correlation = numeric(),
      key     = factor(),
      clusterID   = character(),
      y           = integer()
    )

    plot_df <- data.frame(
      start   = 0,
      y       = 20,
      cluster = factor(0)
    )

  } else {

    corr <- read.table(input_path, header = FALSE, sep = "\t")
    colnames(corr) <- c("chr", "start", "end", "CpG_ID", "correlation", "pval", "score")
    n <- nrow(corr)

    if (n >= 2) {

      pos <- as.numeric(corr$start)
      corr$sign <- ifelse(corr$score > 0, "pos", ifelse(corr$score < 0, "neg", "neut"))
      corr <- corr[order(corr$start), ]
      corr$sign_run <- cumsum(c(TRUE,diff(corr$sign) != 0))
      corr <- corr %>%
	group_by(sign_run) %>%
	mutate(
	       dist = c(0, diff(start)),
	       cluster = cumsum(dist > max_dist) 
	       )
     
      corr$clusterID <- paste0(TxID, "_", corr$sign, "_", corr$sign_run, "_", corr$cluster)
      corr$key <- as.factor(paste0(corr$sign_run, "_", corr$cluster))
      corr$y         <- seq_len(nrow(corr))

      write_df <- corr
      plot_df  <- corr[, c("start", "y","key")]

    } else if (n == 1) {

      corr$key   <- factor(0)
      corr$clusterID <- paste0(TxID, "_0")
      corr$y         <- 1

      write_df <- corr
      plot_df  <- corr[, c("start", "y", "key")]

    } else {
      write_df <- data.frame(
        chr         = character(),
        start       = numeric(),
        end         = numeric(),
        CpG_ID      = character(),
        correlation = numeric(),
        key     = factor(),
        clusterID   = character(),
        y           = integer()
      )

      plot_df <- data.frame(
        start   = 0,
        y       = 20,
        key = factor(0)
      )
    }
  }

  write.table(write_df, file = output_corr, sep = "\t", quote = FALSE, row.names = FALSE)

  pdf(output_plot, width = 10, height = 7)
  print(
    ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = start, y = y, color = key)
    ) +
      ggplot2::geom_point(size = 3) +
      ggplot2::scale_color_brewer(palette = "Set1") +
      ggplot2::labs(
        title = paste0(TxID),
        x = "Genomic Position (bp)",
        y = "",
        color = "Cluster"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.y  = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank()
      )
  )
  dev.off()
}
