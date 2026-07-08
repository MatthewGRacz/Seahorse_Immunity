pca_base <- prcomp(zscores, center = TRUE, scale. = TRUE)

get_dapc_analysis <- function(zscores, num_tests, pipeline_path){
  
  #supertype and cluster used interchangeably here
  
  #report other summary stats of the cluster estimates
  
  cluster_data <- find.clusters(zscores, 
                                n.clust=17,
                                n.pca=ncol(zscores))
  
  #shows BIC graph of clusters, to find optimal number of clusters
  #keep n.pca=ncol(zscores) because you want to use all the variance possible, to overshoot
  #the mode cluster number is used as the number of clusters
  
  test_dapc <- dapc(zscores, 
                    cluster_data$grp, 
                    n.pca=ncol(zscores), 
                    n.da=1000)
  
  #all PCs kept from find.clusters for all amino acids' Z-scores
  #overshooting helps use all values between 1 and the overshot npca number
  #helps to find the most conservative, statistically probable npca number
  #n.da maxes out at n.pca-1, finds the maximum number of axes for the cluster comparisons; 
  #setting n.da too high will cause it to self-correct to the maximum possible value
  
  a_spline_data <- optim.a.score(test_dapc, n.sim=10)
  #finds the optimal number of principal components
  
  final_dapc <- dapc(zscores, 
                     grp = cluster_data$grp, 
                     n.pca = a_spline_data$best,
                     n.da = 1000)
  
  #the groups are the actual separated groups using the best data, accounting for every minor detail
  #DAPC can now get the best analysis of the recombinants using only the optimized number of PCs
  #it looks at the top n.pca separations and keeps those, while the rest are deleted
  
  cat("\n", final_dapc$var * 100, "% of Variance explained using ", final_dapc$n.pca, " PCs and ", length(final_dapc$eig), "/", length(final_dapc$pca.eig), " eigenvalues")
  #how much variance is captured by the optimal number of PCs, 
  #and how many eigenvalues were possible
  
  indv_by_supertype_df <- data.frame(Recombinant = row.names(zscores),
                                     Supertype = final_dapc$assign)
  
  #which cluster is each recombinant in
  

  
  #this gets all of the linear discriminants and plots them against each other in a multi-page plot 
  #helps visualize X-dimensional space, where X is the number of LDs
  
  return(list(final_dapc, a_spline_data))
  
}

num_runs <- 3

results <- replicate(num_runs, get_dapc_analysis(zscores, 1, "Pipeline/"), simplify = FALSE)

View(results)

dapcs <- lapply(1:num_runs, function(i) {
  
  dapc_obj <- results[[i]][[1]] 
  
  df <- data.frame(
    Recombinant = rownames(zscores),
    Assignment = dapc_obj$assign
  )
  
  colnames(df)[2] <- paste0("Run_", i)
  
  return(df)
})

merged_dapc <- Reduce(function(x, y) merge(x, y, by = "Recombinant", suffixes = c("", "_new")), dapcs)
merged_dapc <<- merged_dapc

View(merged_dapc)

ari_scores <- sapply(2:num_runs, function(i) { adjustedRandIndex(merged_dapc[["Run_1"]], merged_dapc[[paste0("Run_", i)]]) })

names(ari_scores) <- paste0("Run_", 2:num_runs)

cat("\n--- True Adjusted Rand Index (ARI) Scores ---\n")

print(ari_scores)
cat("\nMean: ", mean(ari_scores))
cat("\nMedian: ", median(ari_scores))


table(merged_dapc$Run_1, merged_dapc$Run_2)




