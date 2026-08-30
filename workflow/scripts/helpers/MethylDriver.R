
MethylDriver = function(home_region_methyl_change, regions_genes, closest_regions_list, multiple_correction = "fdr", q_value_cutoff = 0.05){
  
  # Check that names home_region_methyl_change vector has names
  if(is.null(names(home_region_methyl_change))){simpleError("home_region_methyl_change should be a named vector")}
  
  # Check that names of home_region_methyl_change match names of closest_regions_list
  if(!all(names(home_region_methyl_change) == names(closest_regions_list))){simpleError("Names of home_region_methyl_change do not match those of closest_regions_list")}
  
  # Check to see that input arguments belong to the set of allowed values
  if(!multiple_correction %in% p.adjust.methods){simpleError("Incorrect value for multiple_correction")}
  
  # Calculate for each region the distribution of the feature of interest for its neighbours
  feature_change_distribution_for_each_region = lapply(closest_regions_list, function(x) sort(home_region_methyl_change[x]))
  
  # Create a data.frame with the results. Add name of regions, genes associated with regions as a list, the value of the feature of interest for the region, 
  # the mean, the standard deviation for the feature of interest for the region's closest neighbours and the mean distance to the regions closest neighbours
  region_normalization_df = data.frame(name = names(home_region_methyl_change), gene = unname(I(regions_genes)), home_region_methyl_change = home_region_methyl_change,
                                       neighborhood_methyl_change_mean = sapply(feature_change_distribution_for_each_region, function(x) mean(x, na.rm = T)),
                                       neighborhood_methyl_change_sd = sapply(feature_change_distribution_for_each_region, sd)) 
  # Add z-score for the region by subtracting the mean of the closest neighbours from the value for a region and dividing by the standard deviation of the neighbours
  region_normalization_df$z_score = (region_normalization_df$home_region_methyl_change - region_normalization_df$neighborhood_methyl_change_mean)/region_normalization_df$neighborhood_methyl_change_sd
  # Test if the z-score differs significantly from the distribution of values for the closest neighbours (either greater than, less than or two-sided)
  region_normalization_df$p_value =  2 * pnorm(-abs(region_normalization_df$z_score))
  
  # Calculate q-value using specified method testing correction method.
  region_normalization_df$q_value = p.adjust(region_normalization_df$p_value, method = multiple_correction)
  
  # Remove closest neighbour distributions from results table 
  region_normalization_df$closest_regions_distribution = NULL
  
  # Make sure region_normalization_df is a dataframe and not a tibble
  region_normalization_df = data.frame(region_normalization_df)
  
  # Arrange results by q_value
  region_normalization_df = dplyr::arrange(region_normalization_df, q_value)
  row.names(region_normalization_df) = NULL
  
  # Get significant hypermethylated and significant hypomethylated genes
  significant_hypermethylated_genes = unique(dplyr::filter(region_normalization_df, q_value < q_value_cutoff, home_region_methyl_change > 0)$gene)
  significant_hypomethylated_genes = unique(dplyr::filter(region_normalization_df, q_value < q_value_cutoff, home_region_methyl_change < 0)$gene)
  
  # Return specified results object
  return(list(hypermethylated_genes = significant_hypermethylated_genes, hypomethylated_genes = significant_hypomethylated_genes, full_results = region_normalization_df, feature_change_distribution_for_each_region = feature_change_distribution_for_each_region))
}