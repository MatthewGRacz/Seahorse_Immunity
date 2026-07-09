library("clue")

setwd("/Users/mattracz/Projects/Wilson_Lab")

get_dapc_analysis <- function(zscores, num_supertypes, pipeline_path){
  
  cluster_data <- find.clusters(zscores, 
                                n.clust=num_supertypes,
                                n.pca=ncol(zscores))
  
  test_dapc <- dapc(zscores, 
                    cluster_data$grp, 
                    n.pca=ncol(zscores), 
                    n.da=1000)
    
  a_spline_data <- optim.a.score(test_dapc, n.sim=1)
  
  final_dapc <- dapc(zscores, 
                     grp = cluster_data$grp, 
                     n.pca = a_spline_data$best,
                     n.da = 1000)
  
  indv_by_supertype_df <- data.frame(Recombinant = row.names(zscores),
                                     Supertype = final_dapc$assign)
  
  return(final_dapc)
  
}
#gets a DAPC for each run, where it assigns supertypes to the recombinants

num_runs <- 5

num_supertypes <- 17

cat("\nRunning DAPC analyses... This may take a moment...")

results <- replicate(num_runs, get_dapc_analysis(zscores, num_supertypes, "Pipeline/"), simplify = FALSE)
#run num_supertypes many DAPC analyses

merged_dapc <- data.frame(
  
  RECOMBINANT = rownames(zscores),
  
  setNames(lapply(results, function(x){ x$assign }), paste0("RUN_", 1:num_runs))
  
)
#all of the DAPC supertypes for each recombinant in a massive dataframe

ari_df2 <- data.frame(
  
  RUN = paste0("RUN_", 2:num_runs),
  
  ARI = sapply(2:num_runs, function(i) { 
    
    adjustedRandIndex(merged_dapc[["RUN_1"]], merged_dapc[[paste0("RUN_", i)]]) 
    
  })
)
#puts into a dataframe the ARI for each run's DAPC compared to run1's 

comparison_cols <- colnames(merged_dapc)[colnames(merged_dapc) != c("RECOMBINANT", "RUN_1")]
#the columns that run1's supertype numbers are compared to, to form the supertype dictionary

translated_cols <- lapply(comparison_cols, function(col_x){
  
  translated_supertypes <- as.numeric(solve_LSAP(table(merged_dapc[["RUN_1"]], merged_dapc[[col_x]]), maximum = TRUE))
  #compares run1 to runX and puts them into a cross-matrix
  #feeds that cross-matrix to the Hungarian Algorithm
  
  supertype_dictionary <- setNames(1:max(as.numeric(levels(merged_dapc[[col_x]]))), translated_supertypes)
  #creates a dictionary based on the Hungarian Algorithm's best supertype numbers match
  
  return(as.numeric(supertype_dictionary[as.character(merged_dapc[[col_x]])]))
  #translates the supertype numbers from runX into run1-equivalent supertypes
  #returns the column of runX's supertypes as run1-translated supertype numbers
  
})

jumper_df <- data.frame(RECOMBINANT=merged_dapc[["RECOMBINANT"]], 
                        INDIVIDUAL=sub("_.*", "", merged_dapc[["RECOMBINANT"]]), 
                        RUN_1=merged_dapc[["RUN_1"]], 
                        setNames(translated_cols, toupper(comparison_cols)))

#a bit hardcoded, but for a better look
#assigns names and makes runX have runX's translated supertypes

stabilities <- lapply(jumper_df$RECOMBINANT, function(recomb){
  
  recomb_data <- as.numeric(jumper_df[jumper_df$RECOMBINANT==recomb, -c(1, 2)])
  
  best_supertype <- as.numeric(names(which.max(table(recomb_data))))
  
  return((sum(recomb_data == best_supertype)-1) * 100 / (length(recomb_data)-1))
  
})

jumper_df$STABILITY <- as.numeric(unlist(stabilities))
#adds up, out of the total number of runs, the number of runs with the same, most frequent supertype number
#the minus 1 excludes run1 comparing with itself 
#this is the optimistic approach, so if run1 was an outlier, and a recombinant is mostly in another supertype, 
#that supertype is deemed the best-match and therefore the percent of the time that the recombinant is within 
#the better supertype is counted instead

recomb_stabilities <- data.frame(jumper_df["RECOMBINANT"], jumper_df["INDIVIDUAL"], jumper_df["STABILITY"])
#summary stats dataframe, with each recombinant, its individual and its stability, next to each other 

indv_stabilities <- aggregate(STABILITY ~ INDIVIDUAL, data = jumper_df, FUN = function(x) {
  c(MEDIAN = median(x, na.rm = TRUE),
    MEAN   = mean(x, na.rm = TRUE),
    RANGE  = diff(range(x, na.rm = TRUE)))
})

indv_stabilities <- cbind(indv_stabilities["INDIVIDUAL"], indv_stabilities$STABILITY)
#gets the median, mean, and range stabilities per individual, for individual-based analyses 

write_csv(jumper_df, paste0(pipeline_path, "jumpers_dataframe.csv"))

write_csv(recomb_stabilities, paste0(pipeline_path, "recombinant_stabilities.csv"))

write_csv(indv_stabilities, paste0(pipeline_path, "individual_stabilities.csv"))

write_csv(ari_df2, paste0(pipeline_path, "ari_df2.csv"))


