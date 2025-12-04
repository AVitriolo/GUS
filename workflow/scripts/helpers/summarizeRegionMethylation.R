summarizeRegionMethylation <- function(meth_rse, assay = 1, genomic_regions, genomic_region_names = NULL, col_summary_function = "colMeans2",
  keep_metadata_cols = FALSE, max_sites_per_chunk = floor(62500000/ncol(meth_rse)), na.rm = TRUE, BPPARAM = BiocParallel::SerialParam(), ...){
  
  # Check that inputs have the correct data type
  stopifnot(is(meth_rse, "RangedSummarizedExperiment"), is(assay, "numeric") | is(assay, "character"),
    is(genomic_regions, "GRanges"), is(genomic_region_names, "character") | is.null(genomic_region_names),
    S4Vectors::isTRUEorFALSE(keep_metadata_cols), 
    (is(max_sites_per_chunk, "numeric") & max_sites_per_chunk >= 1) | is.null(max_sites_per_chunk),
    is(col_summary_function, "character") & length(col_summary_function) == 1, S4Vectors::isTRUEorFALSE(na.rm), 
    is(BPPARAM, "BiocParallelParam"))
  
  # Check that col_summary_function is one of the column summary functions from MatrixGenerics and then extract the function
  if(!col_summary_function %in% ls("package:MatrixGenerics")){
    stop("col_summary_function should be the name of one of the column summary functions from MatrixGenerics e.g. \"colSums2\"")
  }
  col_summary_function = get(col_summary_function, envir = asNamespace("MatrixGenerics"))
  
  # If genomic_region_names is NULL, attempt to use names of genomic_regions
  if(is.null(genomic_region_names)){genomic_region_names <- names(genomic_regions)}
    
  # Add names to genomic_regions if they are not already present and also check that no names are duplicated. 
  if(is.null(genomic_region_names)){
    message("No names for provided regions so using as.character(genomic_regions) as names")
    genomic_region_names <- as.character(genomic_regions)
    names(genomic_regions) <- genomic_region_names
  } else {
    if(length(genomic_region_names) != length(genomic_regions)){
      stop("genomic_region_names must be the same length as genomic_regions")
    } 
    if(anyDuplicated(genomic_region_names)){
      stop("genomic_region_names cannot contain duplicates")
    } else {
      names(genomic_regions) <- genomic_region_names
    }
  }
  
  # Check that there are some seqlevels from genomic_regions present in meth_rse and 
  # give an error if not and a warning if some are missing 
  if(all(!seqlevels(genomic_regions) %in% seqlevels(meth_rse))){
    stop("There are no seqlevels in common between genomic_regions and meth_rse")
  } 
  if(any(!seqlevels(genomic_regions) %in% seqlevels(meth_rse))){
    message("There are some seqlevels from genomic_regions missing from meth_rse")
  } 
  
  # Split genomic regions into chunks based on the number of methylation sites that they cover
  genomic_region_bins <- .chunk_regions(meth_rse = meth_rse, genomic_regions = genomic_regions, 
    max_sites_per_chunk = max_sites_per_chunk, ncores = BiocParallel::bpnworkers(BPPARAM))
  
  # For each sequence get methylation of the associated regions
  region_methylation <- BiocParallel::bpmapply(FUN = .summarize_chunk_methylation, 
    chunk_regions = genomic_region_bins, MoreArgs = list(meth_rse = meth_rse, 
      assay = assay, col_summary_function = col_summary_function, na.rm = na.rm, ...), 
    BPPARAM = BPPARAM, SIMPLIFY = FALSE)

  # Combine data.frames for each chunk
  region_methylation <- dplyr::bind_rows(region_methylation)
  gc()
  
  # Turn rownames into a column and convert the result to a data.table
  region_methylation <- data.table::data.table(tibble::rownames_to_column(region_methylation, "region_name"))
  
  # Create a data.table with the genomic_region_names in the correct order
  genomic_region_names_df <- data.table::data.table(region_name = genomic_region_names)
  
  # Put rows in same order as regions in genomic_regions. Adds rows with NA values for regions which didn't overlap any methylation sites. 
  region_methylation <- data.table::merge.data.table(genomic_region_names_df, region_methylation, 
    by = "region_name", all.x = TRUE, sort = FALSE)

  # Add metadata from genomic_regions if specified
  if(keep_metadata_cols){
    region_methylation <- cbind(region_methylation, data.frame(mcols(genomic_regions)))
  }
  
  # Turn region_name back into row.names
  region_methylation <- data.frame(tibble::column_to_rownames(region_methylation, "region_name"))
  
  # Run the garbage collection and return region_methylation
  gc()
  return(region_methylation)

}