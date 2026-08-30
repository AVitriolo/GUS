# =============================================================================
# make_feature_matrix_from_bigwigs_and_bed()
# get_encode_info() sourced from helpers.
# bigWigAverageOverBed provided on PATH by the conda env.
# =============================================================================
make_feature_matrix_from_bigwigs_and_bed = function(bigwig_dir, genomic_regions_bed,
                                                     bigwig_pattern = NULL,
                                                     exclude_files = character(0),
                                                     covered_bases_only = T,
                                                     num_cores = 4,
                                                     temp_dir = NULL){
  `%do%` = foreach::`%do%`

  # Get paths to bigwig files in specified directory
  bigwig_file_paths = list.files(path = bigwig_dir, pattern = ".bigWig", full.names = T)
  print(bigwig_file_paths)
  # Filter for bigwig files with names matching a specified pattern
  if (!is.null(bigwig_pattern)){
    bigwig_file_paths = grep(bigwig_pattern, bigwig_file_paths, value = T)
  }
  bigwig_file_paths <- bigwig_file_paths[!basename(bigwig_file_paths) %in% exclude_files]
  # Get names of regions from BED file
  region_names = data.table::fread(genomic_regions_bed, header = F, sep = "\t", stringsAsFactors = F)$V4
  print(region_names)
  # Create names for the output files that will be produced by bigWigAverageOverBed
  average_over_bed_file_names = paste0(gsub(".bigWig", "", basename(bigwig_file_paths)), "_average_over_bed.txt")

  # Create a temporary directory to store these output files
  if (is.null(temp_dir)) temp_dir <- tempfile("temp_average_over_bed_results_")
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

  # Use bigWigAverageOverBed to calculate mean value for each bigwig file for each region in the BED file
  parallel::mclapply(
    seq_along(bigwig_file_paths),
    function(bigwig_number) {

      input_bw <- bigwig_file_paths[bigwig_number]

      output_file <- file.path(
        temp_dir,
        average_over_bed_file_names[bigwig_number]
      )

      status <- system2(
        "bigWigAverageOverBed",
        args = c(
          input_bw,
          genomic_regions_bed,
          output_file
        )
      )

      if (!identical(status, 0L)) {
        stop(sprintf(
          "bigWigAverageOverBed failed for %s (exit status %s)",
          input_bw,
          status
        ))
      }

      NULL
    },
    mc.cores = num_cores,
    mc.preschedule = FALSE
  )
  # Select the mean or mean0 columns from the bigWigAverageOverBed output files
  mean_col = ifelse(covered_bases_only, 6, 5)
  result_files <- file.path(temp_dir, average_over_bed_file_names)
  row_counts <- sapply(result_files, function(f) nrow(data.table::fread(f, header=FALSE)))
  print(row_counts)
  # Create a matrix with rows corresponding to genomic regions and columns corresponding to genomic features
  genomic_regions_chromatin_feature_matrix = as.matrix(foreach::foreach(result_file =
                                                                          result_files, .combine = "cbind") %do% {
                                                                            data.table::fread(result_file, sep = "\t", header = F, select = mean_col)
                                                                          })
  # Assign the genomic region names as rownames of the matrix and the feature names as its column names
  rownames(genomic_regions_chromatin_feature_matrix) = region_names
  colnames(genomic_regions_chromatin_feature_matrix) = gsub("_average_over_bed.txt", "", average_over_bed_file_names)
  bw_files <- colnames(genomic_regions_chromatin_feature_matrix)
  # Extract only the accession from the filename
  accessions <- sub("\\.bigWig$", "", bw_files)
  # Query ENCODE for all accessions
  encode_info <- lapply(accessions, get_encode_info)
  targets <- sapply(encode_info, function(x) if (!is.null(x["target"])) x["target"] else NA)
  # Make new column names: "ENCFF003NRZ_H3K27ac"
  new_colnames <- paste(accessions, targets, sep = "_")

  # Assign new column names
  colnames(genomic_regions_chromatin_feature_matrix) <- new_colnames
  return(genomic_regions_chromatin_feature_matrix)
}