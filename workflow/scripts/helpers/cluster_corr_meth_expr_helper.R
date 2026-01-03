cluster_corr <- function(input_path, output_corr, output_plot, TxID, sign) {

  plot_df  <- NULL
  write_df <- NULL

  if (!file.exists(input_path) || file.size(input_path) == 0) {

    write_df <- data.frame(
      chr         = character(),
      start       = numeric(),
      end         = numeric(),
      CpG_ID      = character(),
      correlation = numeric(),
      cluster     = factor(),
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
    colnames(corr) <- c("chr", "start", "end", "CpG_ID", "correlation")
    n <- nrow(corr)

    if (n >= 2) {

      pos <- as.numeric(corr$start)
      dist_mat <- dist(pos, method = "euclidean")
      hc <- hclust(dist_mat, method = "complete")

      clusters_dyn <- dynamicTreeCut::cutreeDynamic(
        dendro = hc,
        distM = as.matrix(dist_mat),
        deepSplit = 2,
        pamRespectsDendro = FALSE,
        minClusterSize = 1
      )

      corr$cluster   <- factor(clusters_dyn)
      corr           <- corr[order(corr$start), ]
      corr$clusterID <- paste0(TxID, "_", corr$cluster, "_", sign)
      corr$y         <- seq_len(nrow(corr))

      write_df <- corr
      plot_df  <- corr[, c("start", "y", "cluster")]

    } else if (n == 1) {

      corr$cluster   <- factor(0)
      corr$clusterID <- paste0(TxID, "_0_", sign)
      corr$y         <- 1

      write_df <- corr
      plot_df  <- corr[, c("start", "y", "cluster")]

    } else {
      write_df <- data.frame(
        chr         = character(),
        start       = numeric(),
        end         = numeric(),
        CpG_ID      = character(),
        correlation = numeric(),
        cluster     = factor(),
        clusterID   = character(),
        y           = integer()
      )

      plot_df <- data.frame(
        start   = 0,
        y       = 20,
        cluster = factor(0)
      )
    }
  }

  write.table(write_df, file = output_corr, sep = "\t", quote = FALSE, row.names = FALSE)

  pdf(output_plot, width = 10, height = 7)
  print(
    ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = start, y = y, color = cluster)
    ) +
      ggplot2::geom_point(size = 3) +
      ggplot2::scale_color_brewer(palette = "Set1") +
      ggplot2::labs(
        title = paste0(TxID, "_", sign),
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
