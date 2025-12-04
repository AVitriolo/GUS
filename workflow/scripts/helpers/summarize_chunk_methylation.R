.summarize_chunk_methylation <- function(chunk_regions, meth_rse, assay, col_summary_function, na.rm, ...){
  
  # Subset meth_rse_for_chunk for regions overlapping chunk_regions
  meth_rse_for_chunk <- subsetByOverlaps(meth_rse, chunk_regions)
  invisible(gc()) 
  
  # Find the overlaps of chunk_regions and meth_rse_for_chunk
  overlaps_df <- data.frame(findOverlaps(chunk_regions, meth_rse_for_chunk))
  
  # Add region names to overlaps_df
  overlaps_df$genomic_region_name <- names(chunk_regions)[overlaps_df$queryHits]
  
  # Create a list matching region names to rows of meth_rse_for_chunk
  region_names_to_rows_list <- split(overlaps_df$subjectHits, overlaps_df$genomic_region_name)
  
  # Read all values from specified assay of meth_rse_for_chunk into memory and run the garbage collection
  meth_values <- as.matrix(SummarizedExperiment::assay(meth_rse_for_chunk, i = assay))

  gc()
  
  # Summarize methylation values 
  meth_summary <- lapply(region_names_to_rows_list, function(x) 
    col_summary_function(meth_values[x, , drop = FALSE], useNames = TRUE, na.rm = na.rm, ...))
  rm(meth_values); gc()
  
  # Combine meth_summary into a single table
  meth_summary <- do.call("rbind", meth_summary)
  
  # Convert meth_summary to a data.frame
  meth_summary <- data.frame(meth_summary)
  gc()
  meth_summary

}
