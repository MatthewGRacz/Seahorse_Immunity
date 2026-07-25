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
    
  a_spline_data <- optim.a.score(test_dapc, n.sim=100)
  
  final_dapc <- dapc(zscores, 
                     grp = cluster_data$grp, 
                     n.pca = a_spline_data$best,
                     n.da = 1000)
  
  indv_by_supertype_df <- data.frame(Recombinant = row.names(zscores),
                                     Supertype = final_dapc$assign)
  
  return(final_dapc)
  
}
#gets a DAPC for each run, where it assigns supertypes to the recombinants

num_runs <- 1200

num_runs <- num_runs + 1
#run1 for a baseline, the other runs as comparisons

num_supertypes <- 17

cat("\nRunning DAPC analyses... This may take a moment...")

results <- replicate(num_runs, get_dapc_analysis(zscores, num_supertypes, "Pipeline/"), simplify = FALSE)
#run num_supertypes many DAPC analyses

merged_dapc <- data.frame(
  
  RECOMBINANT = rownames(zscores),
  
  setNames(lapply(results, function(x){ x$assign }), paste0("RUN_", 1:num_runs))
  
)
#all of the DAPC supertypes for each recombinant in a massive dataframe

ari_df <- data.frame(
  
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

run_cols <- grep("^RUN_", colnames(jumper_df), value = TRUE)
run_cols <- run_cols[run_cols != "RUN_1"]

jumper_df$BEST <- apply(jumper_df[, run_cols], 1, function(recomb_row){
  
  recomb_data <- as.numeric(unlist(recomb_row))
  
  return(as.numeric(names(which.max(table(recomb_data)))))
  
})

jumper_df$STABILITY <- apply(jumper_df[, run_cols], 1, function(recomb_row){
  
  recomb_data <- as.numeric(unlist(recomb_row))
  
  return(100*(max(table(recomb_data))-1)/(length(recomb_data)-1))
  
})

jumper_df$PREFERENCE <- 100 * (jumper_df$STABILITY - (100/num_supertypes)) / (100-(100/num_supertypes))

jumper_df$INSTABILITY <- 100 - jumper_df$STABILITY

#adds up, out of the total number of runs, the number of runs with the same, most frequent supertype number
#the minus 1 excludes run1 comparing with itself 
#this is the optimistic approach, so if run1 was an outlier, and a recombinant is mostly in another supertype, 
#that supertype is deemed the best-match and therefore the percent of the time that the recombinant is within 
#the better supertype is counted instead

jumper_df$UNIQUE_STS <- apply(jumper_df[, run_cols], 1, function(recomb_row){
  
  return(length(unique(as.numeric(recomb_row))))
  
})

jumper_df$EFFECTIVE_STS_HILL_SIMPSON <- apply(jumper_df[, run_cols], 1, function(recomb_row){
  
  freq_table <- sort(table(as.numeric(recomb_row))/length(recomb_row), decreasing = TRUE)
  
  return(1/(sum((freq_table)^2)))
  
})

#Hill Number for species richness using the inverse of the Ginni-Simpson Index
#Since we assume that dominant supertypes are likely where the recombinant is, and that
#less-frequent supertypes are likely noise, this one is used
#The Shannon Hill Number is a bit less exponential, so noise gets weighed more equally

jumper_df$UNEVENNESS <- apply(jumper_df[, run_cols], 1, function(recomb_row){
  
  freq_table <- sort(table(as.numeric(recomb_row))/length(recomb_row), decreasing = TRUE)
  
  if(max(freq_table) == min(freq_table) || length(freq_table)==1) {return(0)}
  
  else{
    
    return(sum(abs(diff(freq_table)))*100)
    
  }
  
})

nonrun_cols <- colnames(jumper_df)[!colnames(jumper_df) %in% run_cols]

jumper_df$NOISE <- jumper_df$UNEVENNESS * jumper_df$UNIQUE_STS / jumper_df$EFFECTIVE_STS_HILL_SIMPSON

jumper_df$NOISE <- jumper_df$NOISE*100/max(jumper_df$NOISE)

#percentage of how much a recombinant travels; 0% is perfect stability, while higher indicates
#that the recombinant has been looped into a greater number of supertypes overall
#different from stability, since a recombinant can have 50% stability among 2 supertypes, with a travelling of 20%,
#or be, among 10 runs, 50% in its most frequent supertype, and 10% among 5 others, hence with a travelling of 50% 

recomb_stabilities <- jumper_df[,nonrun_cols[nonrun_cols != "RUN_1"]]
#summary stats dataframe, with each recombinant, its individual and its stability, next to each other 

indv_stabilities <- aggregate(STABILITY ~ INDIVIDUAL, data = jumper_df, FUN = function(x) {
  c(MEDIAN = median(x, na.rm = TRUE),
    MEAN   = mean(x, na.rm = TRUE),
    RANGE  = diff(range(x, na.rm = TRUE)))
})

indv_stabilities <- cbind(indv_stabilities["INDIVIDUAL"], indv_stabilities)
#gets the median, mean, and range stabilities per individual, for individual-based analyses 

write_csv(jumper_df, paste0(pipeline_path, "jumpers_dataframe.csv"))

write_csv(recomb_stabilities, paste0(pipeline_path, "recombinant_stabilities.csv"))

write_csv(indv_stabilities, paste0(pipeline_path, "individual_stabilities.csv"))

write_csv(ari_df, paste0(pipeline_path, "ari_df.csv"))

GLM_individuals_data <- jumper_df[jumper_df$INDIVIDUAL %in% ordered_names, ]

fixed_individual_supertype_vector <- +(table(1:17 %in% GLM_individuals_data$BEST[indv_x,]))

microbe_attribute_data <- read.csv(paste0(pipeline_path, "GLMOTUSTv2.csv"))
microbe_attribute_data <- column_to_rownames(microbe_attribute_data, "FISH")
JLA_individual_supertype_vector <- microbe_attribute_data[!grepl("^M", colnames(microbe_attribute_data))]
#colname manipulation to add one to each and then pad with 0s the nonexistent ones

JLA_supertypes_col <- lapply(JLA_individual_supertype_vector, function(indv_x){
  
  fixed_individual_supertype_vector <- +(table(1:17 %in% GLM_individuals_data$BEST[indv_x,]))
  #sees if 1:17 is within the best supertypes of this individual
  #repeats for all individuals used in the microbe analysis
  
  translated_supertypes <- as.numeric(solve_LSAP(table(JLA_individual_supertype_vector[indv_x,], fixed_individual_supertype_vector), maximum = TRUE))
  #compares run1 to runX and puts them into a cross-matrix
  #feeds that cross-matrix to the Hungarian Algorithm
  
  supertype_dictionary <- setNames(1:max(as.numeric(levels(fixed_individual_supertype_vector[indv_x,]))), translated_supertypes)
  #creates a dictionary based on the Hungarian Algorithm's best supertype numbers match
  
  return(as.numeric(supertype_dictionary[as.character(fixed_individual_supertype_vector[indv_x,])]))
  #translates the supertype numbers from runX into run1-equivalent supertypes
  #returns the column of runX's supertypes as run1-translated supertype numbers
  
})



