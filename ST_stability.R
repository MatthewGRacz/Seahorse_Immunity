dapcs <- list(read_csv("/Users/mattracz/Projects/Wilson_Lab/Pipeline/dapc1.csv"), 
              read_csv("/Users/mattracz/Projects/Wilson_Lab/Pipeline/dapc2.csv"), 
              read_csv("/Users/mattracz/Projects/Wilson_Lab/Pipeline/dapc3.csv"), 
              read_csv("/Users/mattracz/Projects/Wilson_Lab/Pipeline/dapc4.csv"), 
              read_csv("/Users/mattracz/Projects/Wilson_Lab/Pipeline/dapc5.csv"))
merged_dapc <- Reduce(function(x, y) merge(x, y, by = "RECOMB", suffixes = c("", "_new")), dapcs)

View(merged_dapc)

ari_scores <- c(
  Run2 = adjustedRandIndex(merged_dapc[[2]], merged_dapc[[3]]),
  Run3 = adjustedRandIndex(merged_dapc[[2]], merged_dapc[[4]]),
  Run4 = adjustedRandIndex(merged_dapc[[2]], merged_dapc[[5]]),
  Run5 = adjustedRandIndex(merged_dapc[[2]], merged_dapc[[6]])
)

print(ari_scores)