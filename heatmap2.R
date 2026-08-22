library("Biostrings")
library("seqinr")
library("utils")
library("adegenet")
#remotes::install_version("Peptides", version = 2.4) #version used in analysis
library("Peptides")
library("ggplot2")
library("jsonlite")
library("vegan")
library("tidyverse")
library("rvest")
library("parallel")
library("readxl")
library("mclust")
library("clue")

setwd("/Users/mattracz/Projects/Wilson_Lab")



pipeline_path <- "Pipeline/AB_62426/"

zscores <- read.csv(paste0(pipeline_path, "zscores.csv"))
zscores <- column_to_rownames(zscores, "X")
zscores <<- zscores


get_dapc_analysis <- function(zscores, num_supertypes){
  
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

results <- replicate(num_runs, get_dapc_analysis(zscores, num_supertypes), simplify = FALSE)
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

indv_stabilities <- do.call(data.frame, indv_stabilities)

write_csv(jumper_df, paste0(pipeline_path, "jumpers_dataframe_1200.csv"))

write_csv(recomb_stabilities, paste0(pipeline_path, "recombinant_stabilities_1200.csv"))

write_csv(indv_stabilities, paste0(pipeline_path, "individual_stabilities_1200.csv"))

write_csv(ari_df, paste0(pipeline_path, "ari_df_1200.csv"))



recombs_62426 <- read_recombs("Pipeline/AB_62426/recombs_postRECCO.fasta")
recombs_CP_62426 <- read_recombs("Pipeline/AB_62426/CP/recombs_postRECCO.fasta")



recomb_selection_seqs_62426 <- get_selection_site_seq(recombs_62426, get_datamonkey("Pipeline/AB_62426/slac.csv", 
                                                                                    "Pipeline/AB_62426/meme.json", 
                                                                                    "Pipeline/AB_62426/fel.json"))

CP_recomb_selection_seqs_62426 <- get_selection_site_seq(recombs_CP_62426, get_datamonkey("Pipeline/AB_62426/slac.csv", 
                                                                                          "Pipeline/AB_62426/meme.json", 
                                                                                          "Pipeline/AB_62426/fel.json"))


zscores_62426 <- column_to_rownames(read.csv("Pipeline/AB_62426/zscores.csv"), "X")
jumpers_62426 <- read.csv("Pipeline/AB_62426/jumpers_dataframe_1200.csv")

recomb_jumper_df_62426 <- read.csv("Pipeline/AB_62426/jumpers_dataframe_1200.csv")

microbe_attribute_data <- read.csv(paste0("/Users/mattracz/Projects/Wilson_Lab/Pipeline/GLMOTUSTv2.csv"))



ST_AB4 <- ST_AB[,c("RECOMBINANT", "UNIQUE_RECOMBINANT", "INDIVIDUAL", "BEST", "STABILITY", "JLA")]
ST_AB4$SEQ <- recomb_selection_seqs_62426[match(ST_AB4$RECOMBINANT, names(recomb_selection_seqs_62426))]
ST_AB_CP_62426 <- read.csv("/Users/mattracz/Projects/Wilson_Lab/Pipeline/AB_62426/CP/recombinant_stabilities_1200.csv")
ST_AB_CP_62426 <- ST_AB_CP_62426[,c("RECOMBINANT", "INDIVIDUAL", "BEST")]
ST_AB_CP_62426$CP_SEQ <- CP_recomb_selection_seqs_62426[match(ST_AB_CP_62426$RECOMBINANT, names(CP_recomb_selection_seqs_62426))]
ST_AB4$CP_BEST <- NA
ST_AB4$CP_BEST <- ST_AB_CP_62426$BEST[match(ST_AB4$SEQ, ST_AB_CP_62426$CP_SEQ)]

ST_AB4$PROTEIN <- get_translation(recomb_selection_seqs_62426)[match(ST_AB4$RECOMBINANT, names(recomb_selection_seqs_62426))]
ST_AB4$PROTEIN <- gsub("[^A-Z]", "", ST_AB4$PROTEIN)
ST_AB_CP_62426$PROTEIN <- get_translation(CP_recomb_selection_seqs_62426)[match(ST_AB_CP_62426$RECOMBINANT, names(CP_recomb_selection_seqs_62426))]
ST_AB_CP_62426$PROTEIN <- gsub("[^A-Z]", "", ST_AB_CP_62426$PROTEIN)

ST_AB4$CP_BEST <- ST_AB_CP_62426$BEST[match(ST_AB4$PROTEIN, ST_AB_CP_62426$PROTEIN)]

View(ST_AB_CP_62426)
View(ST_AB4)





parent_supertypes_62426 <- as.factor(jumpers_62426$BEST)

parent_DAPC_62426 <- dapc(zscores_62426, 
                          grp = parent_supertypes_62426, 
                          n.pca = ncol(zscores_62426), 
                          n.da = length(levels(parent_supertypes_62426))-1)


prediction_62426 <- predict(parent_DAPC_62426, 
                            newdata = get_z_scores(strsplit(ST_AB_CP_62426$PROTEIN, "")))

ST_AB_CP_62426$TRANSLATED_BEST <- sprintf("S%02d", as.numeric(as.character(prediction_62426$assign)) - 1)

View(ST_AB_CP_62426)

recombs_data_62426 <- recomb_jumper_df_62426[recomb_jumper_df_62426$INDIVIDUAL %in% ordered_names, ]

#have some line which filters based on noise, like 
#recombs_data_62426 <- recomb_jumper_df_62426[recomb_jumper_df_62426$NOISE < 95] 
#and play with number then generate matrix, heatmap; see which relations hold/dont


matrix_62426 <- table(recombs_data_62426$INDIVIDUAL, factor(recombs_data_62426$BEST, levels = 1:17))
colnames(matrix_62426) <- sprintf("S%02d", 0:16)
matrix_62426 <- matrix_62426[ordered_names,]
matrix_62426 <- +(matrix_62426 > 0)
matrix_62426 <- unclass(matrix_62426)
class(matrix_62426) <- "matrix"


microbe_attribute_data <- column_to_rownames(microbe_attribute_data, "FISH")
JLA_matrix <- microbe_attribute_data[!grepl("^M", colnames(microbe_attribute_data))]
microbe_attribute_data <- microbe_attribute_data[!grepl("^S", colnames(microbe_attribute_data))]
JLA_matrix <- as.matrix(JLA_matrix)

JLA_microbes <- c('M0008', 'M0003', 'M0006', 'M0016', 'M0044', 'M0018', 'M0061', 'M0081', 'M0031', 'M0014', 'M0012', 'M0007', 'M0015', 'M0019', 'M0011', 'M0005', 'M0010', 'M0030', 'M0041', 'M0085', 'M0029', 'M0020', 'M0048', 'M0013', 'M0002', 'M0004', 'M0023', 'M0032', 'M0064', 'M0001', 'M0017', 'M0009')
#also can just use the top ones with 100+ reads
#JLA_microbes <- colnames(microbe_attribute_data[colSums(microbe_attribute_data) > 99])
#play and see what happens
#can also get top 32 of those instead, so
#JLA_microbes <- colnames(sort(microbe_attribute_data[colSums(microbe_attribute_data) > 99], descending=TRUE))[1:32]


microbe_attribute_data_62426 <- cbind(microbe_attribute_data, as.data.frame(matrix_62426))



get_glm_and_heatmap <- function(pipeline_path, 
                                microbe_attribute_data, 
                                JLA_microbes,
                                attribute_symbol, 
                                attribute_name, 
                                attribute_data, 
                                attribute_freqs, 
                                heatmap_title, 
                                heatmap_file_name, 
                                x_scale,
                                y_scale){
  
  attribute_symbol <- paste0("^", attribute_symbol)
  #for regular expression (regex) functions
  
  microbe_data <- microbe_attribute_data[, colnames(microbe_attribute_data) %in% JLA_microbes]
  #microbe data is the same for the G and H (Microbe vs Supertype and Microbe vs Unique_Recombinant)
  
  attribute_data <- attribute_data[rownames(microbe_data), , drop = FALSE]
  #puts it in the same order as JLA
  
  if(is.null(attribute_data)){
    
    #set as NULL for JLA, since her data contains the logical vectors of presence/absence 
    #of individuals' recombinants in supertypes or unique_recombs or whatever other attribute
    #we can only know the abundance of my attribute, since JLA's frequencies are not available
    #and her logical vector can make it so an individual with 4 recombs all of one attribute 
    #get condensed to a 1 for that attribute, thereby losing data
    
    attribute_data <- microbe_attribute_data[, grep(x = colnames(microbe_attribute_data), pattern = attribute_symbol, value=TRUE)]
    
    #this gets the columns of logical vectors for presence/absence for that attribute by individual
    
  }
  
  microbe_attribute_combos <- setNames(
    expand.grid(colnames(microbe_data), 
                colnames(attribute_data), 
                stringsAsFactors = FALSE),
    
    c("Microbe", attribute_name))
  
  #creates every Microbe-Attribute combo to run a GLM on
  
  cat("\nCalcuating GLM... This may take a moment...\n")
  
  cl <- makeCluster(detectCores() - 1)
  #uses all except 1 computer CPU core to run, since can take up to 8 hours!
  
  clusterExport(cl, c("attribute_data", "microbe_data", "attribute_name"), envir = environment())
  #exports relevant data to each subordinate R session that runs in parallel (parApply next line),
  #each of which creates its own mini-environment
  
  GLMresults <- suppressWarnings(parApply(cl, microbe_attribute_combos, MARGIN = 1, FUN = function(microbe_attribute) {
    
    attribute <- unname(microbe_attribute[[attribute_name]])
    microbe   <- unname(microbe_attribute["Microbe"])
    
    #glm_analysis <- suppressWarnings(coef(summary(glm(attribute_data[, attribute] ~ microbe_data[, microbe], family = quasipoisson(link = "log")))))
    
    glm_analysis <- suppressWarnings(coef(summary(glm(attribute_data[, attribute] ~ microbe_data[, microbe], family = quasibinomial(link = "logit")))))
    
    
    #the actual GLM analysis for each association
    
    data.frame(
      RECOMB_MICROBE = paste0(attribute, microbe), 
      MICROBE = microbe, 
      ATTRIBUTE = attribute, 
      SLOPE = glm_analysis[2, 1],
      SE = glm_analysis[2, 2],
      WALDZ = glm_analysis[2, 1] / glm_analysis[2, 2],
      P = glm_analysis[2, 4],
      stringsAsFactors = FALSE
    )
    
  }))
  
  stopCluster(cl)
  #return to normal CPU allocation and function after the parallel computations are done
  
  GLMresults <- do.call(rbind, GLMresults)
  
  GLMresults$SEQ_BONFERRONI_P <- p.adjust(GLMresults$P, method = "holm")
  #Sequential Bonferroni correction (Holm-Bonferroni) along the entire dataframe
  
  GLMresults$SIGNIFICANCE <- as.character(symnum(GLMresults$SEQ_BONFERRONI_P, 
                                                 corr = FALSE, 
                                                 na = FALSE, 
                                                 cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                                 symbols = c("***", "**", "*", "")))
  
  #the corrected p value per GLM analysis gets marked with stars if p <= certain thresholds
  
  GLMresults$MICROBE <- factor(GLMresults$MICROBE, levels = JLA_microbes)
  #keeps the microbes in the order they're in
  
  cat("\nGLM calculations complete!\n")
  
  min_z <- min(GLMresults$WALDZ, na.rm = TRUE)
  max_z <- max(GLMresults$WALDZ, na.rm = TRUE)
  mid_z <- (min_z + max_z)/2
  
  heatmap <- suppressWarnings(ggplot(GLMresults, aes(x = paste0(ATTRIBUTE, "\n",  attribute_freqs[ATTRIBUTE]), y = MICROBE, fill = WALDZ)) +
                                geom_tile(color = "black", size=0.4) +
                                geom_text(aes(label = SIGNIFICANCE, color = WALDZ > mid_z), size = 4, vjust = 0.7)  +
                                scale_color_manual(values = c("TRUE" = "white", "FALSE" = "white"), guide = "none") +
                                scale_fill_gradient2(
                                  #low = "white", mid="#A6A6A6", high = "black", 
                                  low = "blue", mid="hotpink", high = "#ff0026", 
                                  midpoint=mid_z,
                                  limits = c(min_z, max_z),
                                  breaks = c(min_z, max_z),
                                  labels = ceiling(c(min_z, max_z)), 
                                  name = "Wald's Z"
                                ) +
                                guides(fill = guide_colorbar(frame.color = "black", frame.linewidth = 0.2)) +
                                theme_minimal() +
                                scale_x_discrete(position = "top") +
                                theme(
                                  axis.text.x = element_text(size = x_scale * 22), 
                                  axis.text.y = element_text(size = y_scale * 30), 
                                  plot.title = element_text(size = 60 * x_scale, face = "bold"), 
                                  axis.title.x = element_text(size = x_scale * 50, face = "bold"),
                                  axis.title.y = element_text(size = y_scale * 50, face = "bold")
                                ) +
                                labs(title = heatmap_title, x = attribute_name , y = "Microbe"))
  
  #suppressWarnings(ggsave(paste0(pipeline_path, heatmap_file_name, ".pdf"), plot = heatmap, width = x_scale * length(unique(GLMresults$ATTRIBUTE)), height = num_microbes * y_scale, limitsize = FALSE))
  
  print(heatmap)
  
  return(GLMresults)
  
}

pipeline_path <- "Pipeline/AB_62426/"


JLA_MST_data_62426 <- get_glm_and_heatmap(pipeline_path, 
                                          microbe_attribute_data_62426, 
                                          JLA_microbes,
                                          "S", 
                                          "Supertype", 
                                          NULL, 
                                          NULL, 
                                          "JLA Microbe–Supertype Associations 6/24/2026", 
                                          "JLA_Supertype_Microbe_WaldZ_Heatmap",
                                          0.3,
                                          0.3)








library("Biostrings")
library("seqinr")
library("utils")
library("adegenet")
#remotes::install_version("Peptides", version = 2.4) #version used in analysis
library("Peptides")
library("ggplot2")
library("jsonlite")
library("vegan")
library("tidyverse")
library("rvest")
library("parallel")
library("readxl")
library("mclust")
library("clue")

setwd("/Users/mattracz/Projects/Wilson_Lab")



pipeline_path <- "Pipeline/AB_72726/"

zscores <- read.csv(paste0(pipeline_path, "zscores.csv"))
zscores <- column_to_rownames(zscores, "X")
zscores <<- zscores


get_dapc_analysis <- function(zscores, num_supertypes){
  
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

results <- replicate(num_runs, get_dapc_analysis(zscores, num_supertypes), simplify = FALSE)
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

indv_stabilities <- do.call(data.frame, indv_stabilities)

write_csv(jumper_df, paste0(pipeline_path, "jumpers_dataframe_1200.csv"))

write_csv(recomb_stabilities, paste0(pipeline_path, "recombinant_stabilities_1200.csv"))

write_csv(indv_stabilities, paste0(pipeline_path, "individual_stabilities_1200.csv"))

write_csv(ari_df, paste0(pipeline_path, "ari_df_1200.csv"))





recombs_72726 <- read_recombs("Pipeline/AB_72726/recombs_postRECCO.fasta")
recombs_CP_72726 <- read_recombs("Pipeline/AB_72726/CP/recombs_postRECCO.fasta")


recomb_selection_seqs_72726 <- get_selection_site_seq(recombs_72726, get_datamonkey("Pipeline/AB_72726/slac.csv", 
                                                                                    "Pipeline/AB_72726/meme.json", 
                                                                                    "Pipeline/AB_72726/fel.json"))

CP_recomb_selection_seqs_72726 <- get_selection_site_seq(recombs_CP_72726, get_datamonkey("Pipeline/AB_72726/slac.csv", 
                                                                                          "Pipeline/AB_72726/meme.json", 
                                                                                          "Pipeline/AB_72726/fel.json"))

zscores_72726 <- column_to_rownames(read.csv("Pipeline/AB_72726/zscores.csv"), "X")
jumpers_72726 <- read.csv("Pipeline/AB_72726/jumpers_dataframe_1200.csv")

recomb_jumper_df_72726 <- read.csv("Pipeline/AB_72726/jumpers_dataframe_1200.csv")

microbe_attribute_data <- read.csv(paste0("/Users/mattracz/Projects/Wilson_Lab/Pipeline/GLMOTUSTv2.csv"))



ST_AB5 <- ST_AB[,c("RECOMBINANT", "UNIQUE_RECOMBINANT", "INDIVIDUAL", "BEST", "STABILITY", "JLA")]
ST_AB5$SEQ <- recomb_selection_seqs_72726[match(ST_AB5$RECOMBINANT, names(recomb_selection_seqs_72726))]
ST_AB_CP_72726 <- read.csv("/Users/mattracz/Projects/Wilson_Lab/Pipeline/AB_72726/CP/recombinant_stabilities_1200.csv")
ST_AB_CP_72726 <- ST_AB_CP_72726[,c("RECOMBINANT", "INDIVIDUAL", "BEST")]
ST_AB_CP_72726$CP_SEQ <- CP_recomb_selection_seqs_72726[match(ST_AB_CP_72726$RECOMBINANT, names(CP_recomb_selection_seqs_72726))]
ST_AB5$CP_BEST <- NA
ST_AB5$CP_BEST <- ST_AB_CP_72726$BEST[match(ST_AB5$SEQ, ST_AB_CP_72726$CP_SEQ)]

ST_AB5$PROTEIN <- get_translation(recomb_selection_seqs_72726)[match(ST_AB5$RECOMBINANT, names(recomb_selection_seqs_72726))]
ST_AB5$PROTEIN <- gsub("[^A-Z]", "", ST_AB5$PROTEIN)
ST_AB_CP_72726$PROTEIN <- get_translation(CP_recomb_selection_seqs_72726)[match(ST_AB_CP_72726$RECOMBINANT, names(CP_recomb_selection_seqs_72726))]
ST_AB_CP_72726$PROTEIN <- gsub("[^A-Z]", "", ST_AB_CP_72726$PROTEIN)

ST_AB5$CP_BEST <- ST_AB_CP_72726$BEST[match(ST_AB5$PROTEIN, ST_AB_CP_72726$PROTEIN)]

View(ST_AB_CP_72726)
View(ST_AB5)





parent_supertypes_72726 <- as.factor(jumpers_72726$BEST)

parent_DAPC_72726 <- dapc(zscores_72726, 
                          grp = parent_supertypes_72726, 
                          n.pca = ncol(zscores_72726), 
                          n.da = length(levels(parent_supertypes_72726))-1)


prediction_72726 <- predict(parent_DAPC_72726, 
                            newdata = get_z_scores(strsplit(ST_AB_CP_72726$PROTEIN, "")))

ST_AB_CP_72726$TRANSLATED_BEST <- sprintf("S%02d", as.numeric(as.character(prediction_72726$assign)) - 1)

View(ST_AB_CP_72726)

recombs_data_72726 <- recomb_jumper_df_72726[recomb_jumper_df_72726$INDIVIDUAL %in% ordered_names, ]
matrix_72726 <- table(recombs_data_72726$INDIVIDUAL, factor(recombs_data_72726$BEST, levels = 1:17))
colnames(matrix_72726) <- sprintf("S%02d", 0:16)
matrix_72726 <- matrix_72726[ordered_names,]
matrix_72726 <- +(matrix_72726 > 0)
matrix_72726 <- unclass(matrix_72726)
class(matrix_72726) <- "matrix"


microbe_attribute_data <- column_to_rownames(microbe_attribute_data, "FISH")
JLA_matrix <- microbe_attribute_data[!grepl("^M", colnames(microbe_attribute_data))]
microbe_attribute_data <- microbe_attribute_data[!grepl("^S", colnames(microbe_attribute_data))]
JLA_matrix <- as.matrix(JLA_matrix)


microbe_attribute_data_72726 <- cbind(microbe_attribute_data, as.data.frame(matrix_72726))



get_glm_and_heatmap <- function(pipeline_path, 
                                microbe_attribute_data, 
                                JLA_microbes,
                                attribute_symbol, 
                                attribute_name, 
                                attribute_data, 
                                attribute_freqs, 
                                heatmap_title, 
                                heatmap_file_name, 
                                x_scale,
                                y_scale){
  
  attribute_symbol <- paste0("^", attribute_symbol)
  #for regular expression (regex) functions
  
  microbe_data <- microbe_attribute_data[, colnames(microbe_attribute_data) %in% JLA_microbes]
  #microbe data is the same for the G and H (Microbe vs Supertype and Microbe vs Unique_Recombinant)
  
  attribute_data <- attribute_data[rownames(microbe_data), , drop = FALSE]
  #puts it in the same order as JLA
  
  if(is.null(attribute_data)){
    
    #set as NULL for JLA, since her data contains the logical vectors of presence/absence 
    #of individuals' recombinants in supertypes or unique_recombs or whatever other attribute
    #we can only know the abundance of my attribute, since JLA's frequencies are not available
    #and her logical vector can make it so an individual with 4 recombs all of one attribute 
    #get condensed to a 1 for that attribute, thereby losing data
    
    attribute_data <- microbe_attribute_data[, grep(x = colnames(microbe_attribute_data), pattern = attribute_symbol, value=TRUE)]
    
    #this gets the columns of logical vectors for presence/absence for that attribute by individual
    
  }
  
  microbe_attribute_combos <- setNames(
    expand.grid(colnames(microbe_data), 
                colnames(attribute_data), 
                stringsAsFactors = FALSE),
    
    c("Microbe", attribute_name))
  
  #creates every Microbe-Attribute combo to run a GLM on
  
  cat("\nCalcuating GLM... This may take a moment...\n")
  
  cl <- makeCluster(detectCores() - 1)
  #uses all except 1 computer CPU core to run, since can take up to 8 hours!
  
  clusterExport(cl, c("attribute_data", "microbe_data", "attribute_name"), envir = environment())
  #exports relevant data to each subordinate R session that runs in parallel (parApply next line),
  #each of which creates its own mini-environment
  
  GLMresults <- suppressWarnings(parApply(cl, microbe_attribute_combos, MARGIN = 1, FUN = function(microbe_attribute) {
    
    attribute <- unname(microbe_attribute[[attribute_name]])
    microbe   <- unname(microbe_attribute["Microbe"])
    
    #glm_analysis <- suppressWarnings(coef(summary(glm(attribute_data[, attribute] ~ microbe_data[, microbe], family = quasipoisson(link = "log")))))
    
    glm_analysis <- suppressWarnings(coef(summary(glm(attribute_data[, attribute] ~ microbe_data[, microbe], family = quasibinomial(link = "logit")))))
    
    
    #the actual GLM analysis for each association
    
    data.frame(
      RECOMB_MICROBE = paste0(attribute, microbe), 
      MICROBE = microbe, 
      ATTRIBUTE = attribute, 
      SLOPE = glm_analysis[2, 1],
      SE = glm_analysis[2, 2],
      WALDZ = glm_analysis[2, 1] / glm_analysis[2, 2],
      P = glm_analysis[2, 4],
      stringsAsFactors = FALSE
    )
    
  }))
  
  stopCluster(cl)
  #return to normal CPU allocation and function after the parallel computations are done
  
  GLMresults <- do.call(rbind, GLMresults)
  
  GLMresults$SEQ_BONFERRONI_P <- p.adjust(GLMresults$P, method = "holm")
  #Sequential Bonferroni correction (Holm-Bonferroni) along the entire dataframe
  
  GLMresults$SIGNIFICANCE <- as.character(symnum(GLMresults$SEQ_BONFERRONI_P, 
                                                 corr = FALSE, 
                                                 na = FALSE, 
                                                 cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                                 symbols = c("***", "**", "*", "")))
  
  #the corrected p value per GLM analysis gets marked with stars if p <= certain thresholds
  
  GLMresults$MICROBE <- factor(GLMresults$MICROBE, levels = JLA_microbes)
  #keeps the microbes in the order they're in
  
  cat("\nGLM calculations complete!\n")
  
  min_z <- min(GLMresults$WALDZ, na.rm = TRUE)
  max_z <- max(GLMresults$WALDZ, na.rm = TRUE)
  mid_z <- (min_z + max_z)/2
  
  heatmap <- suppressWarnings(ggplot(GLMresults, aes(x = paste0(ATTRIBUTE, "\n",  attribute_freqs[ATTRIBUTE]), y = MICROBE, fill = WALDZ)) +
                                geom_tile(color = "black", size=0.4) +
                                geom_text(aes(label = SIGNIFICANCE, color = WALDZ > mid_z), size = 4, vjust = 0.7)  +
                                scale_color_manual(values = c("TRUE" = "white", "FALSE" = "white"), guide = "none") +
                                scale_fill_gradient2(
                                  #low = "white", mid="#A6A6A6", high = "black", 
                                  low = "blue", mid="hotpink", high = "#ff0026", 
                                  midpoint=mid_z,
                                  limits = c(min_z, max_z),
                                  breaks = c(min_z, max_z),
                                  labels = ceiling(c(min_z, max_z)), 
                                  name = "Wald's Z"
                                ) +
                                guides(fill = guide_colorbar(frame.color = "black", frame.linewidth = 0.2)) +
                                theme_minimal() +
                                scale_x_discrete(position = "top") +
                                theme(
                                  axis.text.x = element_text(size = x_scale * 22), 
                                  axis.text.y = element_text(size = y_scale * 30), 
                                  plot.title = element_text(size = 60 * x_scale, face = "bold"), 
                                  axis.title.x = element_text(size = x_scale * 50, face = "bold"),
                                  axis.title.y = element_text(size = y_scale * 50, face = "bold")
                                ) +
                                labs(title = heatmap_title, x = attribute_name , y = "Microbe"))
  
  #suppressWarnings(ggsave(paste0(pipeline_path, heatmap_file_name, ".pdf"), plot = heatmap, width = x_scale * length(unique(GLMresults$ATTRIBUTE)), height = num_microbes * y_scale, limitsize = FALSE))
  
  print(heatmap)
  
  return(GLMresults)
  
}



JLA_MST_data_72726 <- get_glm_and_heatmap(pipeline_path, 
                                          microbe_attribute_data_72726, 
                                          JLA_microbes,
                                          "S", 
                                          "Supertype", 
                                          NULL, 
                                          NULL, 
                                          "JLA Microbe–Supertype Associations 7/27/2026", 
                                          "JLA_Supertype_Microbe_WaldZ_Heatmap",
                                          0.3,
                                          0.3)






