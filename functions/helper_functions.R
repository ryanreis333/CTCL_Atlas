## Get Variable Features for SCT Assay List
get_top_sct_features_across_samples <- function(
    sct_model_list,
    n_top_per_sample,
    m_top_overall,
    residual_variance_col_name = "residual_variance") { # Removed feature_col_name arg
  
  if (!is.list(sct_model_list) || length(sct_model_list) == 0) {
    stop("sct_model_list must be a non-empty list of SCTModel objects.")
  }
  
  if (n_top_per_sample <= 0 || !is.numeric(n_top_per_sample)) {
    stop("n_top_per_sample must be a positive integer.")
  }
  
  if (m_top_overall <= 0 || !is.numeric(m_top_overall)) {
    stop("m_top_overall must be a positive integer.")
  }
  
  all_top_features <- list()
  
  # Step 1: Order every SCTModel layer by residual variance and pull the top 'n' from each
  for (sample_name in names(sct_model_list)) {
    sct_model <- sct_model_list[[sample_name]]
    
    # Ensure it's an SCTModel (or at least has the expected structure)
    if (!inherits(sct_model, "SCTModel") && !all(c("feature.attributes") %in% names(sct_model))) {
      warning(paste("Skipping", sample_name, "as it does not appear to be a valid SCTModel object or lacks @feature.attributes."))
      next
    }
    
    feature_attributes <- sct_model@feature.attributes
    
    # Check if residual_variance column exists
    if (!residual_variance_col_name %in% colnames(feature_attributes)) {
      warning(paste0("Column '", residual_variance_col_name, "' not found in feature.attributes for sample '", sample_name, "'. Skipping."))
      next
    }
    
    # Convert rownames to a column named 'feature'
    feature_attributes_df <- feature_attributes %>%
      rownames_to_column("feature")
    
    # Order by residual variance in descending order
    # Add a rank column for later use
    ranked_features <- feature_attributes_df %>% # Use the df with rownames as a column
      arrange(desc(.data[[residual_variance_col_name]])) %>%
      mutate(rank = row_number()) %>%
      select(feature, rank) # Now 'feature' column exists
    
    # Pull the top 'n_top_per_sample' features
    top_n_features <- ranked_features %>%
      slice_head(n = n_top_per_sample) %>%
      mutate(sample = sample_name) # Add sample identifier
    
    all_top_features[[sample_name]] <- top_n_features
  }
  
  if (length(all_top_features) == 0) {
    message("No valid SCTModel layers were processed. Returning empty data frame.")
    return(data.frame(
      feature = character(0),
      overall_rank = numeric(0),
      frequency = numeric(0),
      avg_individual_rank = numeric(0)
    ))
  }
  
  # Combine all top features from all samples
  combined_top_features <- bind_rows(all_top_features)
  
  # Step 2: Find the top 'm' features that are highest ranked across all samples
  # Calculate frequency and average rank for each feature
  overall_feature_ranking <- combined_top_features %>%
    group_by(feature) %>% # 'feature' is now a consistent column name
    summarise(
      frequency = n(),
      avg_individual_rank = mean(rank),
      .groups = "drop"
    ) %>%
    # Rank features based on frequency (higher is better), then average rank (lower is better)
    arrange(desc(frequency), avg_individual_rank) %>%
    mutate(overall_rank = row_number()) %>%
    slice_head(n = m_top_overall)
  
  return(overall_feature_ranking)
}









