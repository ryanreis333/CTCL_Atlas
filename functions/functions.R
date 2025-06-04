# Function to compute clustering and summary statistics
cluster_summary <- function(mat, k_values) {
  results <- data.frame(
    k = integer(),
    withinss = numeric(),
    silhouette = numeric()
  )
  
  for (k in k_values) {
    # Perform k-means clustering
    kmeans_result <- kmeans(mat, centers = k, nstart = 25, iter.max = 25)
    
    # Compute within-cluster sum of squares
    withinss <- sum(kmeans_result$withinss)
    
    # Compute silhouette score
    dist_matrix <- dist(mat)
    sil <- silhouette(kmeans_result$cluster, dist_matrix)
    silhouette_avg <- mean(sil[, "sil_width"])
    
    # Append results
    results <- rbind(results, data.frame(
      k = k,
      withinss = withinss,
      silhouette = silhouette_avg
    ))
  }
  
  return(results)
}

# Function to extract the top 50 gene names (columns) per cluster based on summed expression values
get_top_50_per_cluster <- function(mat_with_clusters, cluster_id) {
  # Subset rows (cells) belonging to the specified cluster
  cluster_data <- mat_with_clusters[mat_with_clusters$cluster == cluster_id, ]
  
  # Remove the 'cluster' column so that only gene expression data remain
  gene_data <- cluster_data[, setdiff(names(cluster_data), "cluster")]
  
  # Calculate the sum of expression for each gene (i.e., column sums over all cells in the cluster)
  gene_sums <- colSums(gene_data)
  
  # Sort the genes by summed expression in descending order and take the top 50 gene names
  top_50 <- names(sort(gene_sums, decreasing = TRUE)[1:25])
  
  return(top_50)
}