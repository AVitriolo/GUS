library(RANN)

most_similar_n_regions <- function(feature_table, ranges = NULL, n = 100, n_pcs = 20) {
  
  # PCA
  pca_result <- prcomp(
    as.matrix(feature_table),
    scale. = TRUE,
    center = TRUE
  )$x[, seq_len(min(n_pcs, ncol(feature_table)))]
  
  message("donePCA")
  
  # kNN 
  nn <- nn2(
    data = pca_result,
    query = pca_result,
    k = n + 1  # +1 because the closest point is itself
  )
  message("doneknn")
  
  idx <- nn$nn.idx[, -1]  # drop self-match
  
  closest_regions_list <- setNames(
    lapply(seq_len(nrow(idx)), function(i) idx[i, ]),
    rownames(pca_result)
  )
  
  # remove overlapping regions
  if (!is.null(ranges)) {
    overlaps <- GenomicRanges::findOverlaps(
      ranges, ranges, ignore.strand = TRUE
    )
    
    overlap_list <- split(
      S4Vectors::subjectHits(overlaps),
      S4Vectors::queryHits(overlaps)
    )
    
    closest_regions_list <- lapply(seq_along(closest_regions_list), function(i) {
      setdiff(closest_regions_list[[i]], overlap_list[[as.character(i)]])
    })
  }
  
  closest_regions_list
}