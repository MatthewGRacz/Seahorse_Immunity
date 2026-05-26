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

setwd("/Users/mattracz/Projects/Wilson_Lab")



get_recombs <- function(alpha_phasepath, beta_phasepath){
  
  alpha_phasepath <- paste0(alpha_phasepath, "phased.fasta")
  beta_phasepath <- paste0(beta_phasepath, "phased.fasta")
  #file with the phased alleles
  
  alpha_seqs <- read.fasta(alpha_phasepath, as.string = TRUE)
  beta_seqs <- read.fasta(beta_phasepath, as.string = TRUE)
  
  if(sum(names(alpha_seqs) != names(beta_seqs)) == 0){cat("The alpha and beta alleles are a 1:1 match! The pipeline is good to go!\n\n")}
  #all sequence names are the same 
  
  indv_names <- unique(gsub("[ab]$", "", names(alpha_seqs)))
  #get names of individuals whose sequences you're looking at
  
  a1b1 <- toupper(paste0(alpha_seqs[paste0(indv_names, "a")], substr(beta_seqs[paste0(indv_names, "a")], 3, 272)))
  #a = 1st allele, b = 2nd allele, per individual
  #the snipped betas have 2 extra bps ahead and 1 extra bp behind, so cut them out
  #trims to functional reading frames
  a1b2 <- toupper(paste0(alpha_seqs[paste0(indv_names, "a")], substr(beta_seqs[paste0(indv_names, "b")], 3, 272)))
  a2b1 <- toupper(paste0(alpha_seqs[paste0(indv_names, "b")], substr(beta_seqs[paste0(indv_names, "a")], 3, 272)))
  a2b2 <- toupper(paste0(alpha_seqs[paste0(indv_names, "b")], substr(beta_seqs[paste0(indv_names, "b")], 3, 272)))
  
  return(c(
    
    setNames(a1b1, paste0(indv_names, "_a1b1")),
    setNames(a1b2, paste0(indv_names, "_a1b2")),
    setNames(a2b1, paste0(indv_names, "_a2b1")),
    setNames(a2b2, paste0(indv_names, "_a2b2"))
    
  ))
  
}

remove_recco <- function(recombs, pipeline_path){
  
  recco_path <- paste0(pipeline_path, "recco_results.csv")
  
  recco_names <- read.delim(recco_path)$Sequence
  #RECCO's output file is TSV, not CSV like it says
  
  return(recombs[!names(recombs) %in% recco_names])
  
}

get_translation <- function(recombs){
  
  recomb_proteins <- lapply(recombs, function(x) {translate(s2c(x))})
  
  return(recomb_proteins[!sapply(recomb_proteins, function(x) "*" %in% x)])
  #if any have a stop codon by chance of such bps being concatenated, filter them out
  
}

get_datamonkey <- function(slac, meme, fel){
  
  slac <- read.csv(slac)
  
  slac_sites <- c(slac[slac$"P..dN.dS...1." < 0.05, ]$"Site")
  #CSV compresses headers; this corresponds with the positive selection header
  #P [dN/dS > 1]
  
  meme <- fromJSON(meme)$MLE$content$`0`
  #where the data lies
  
  meme_sites <- c(which(meme[,7] < 0.05))
  #7th column is p values, and all p values < 0.05 are pos selection
  
  fel <- fromJSON(fel)$MLE$content$`0`
  
  fel_sites <- c(which(fel[,5] < 0.05))
  #5th column is p values, and all p values < 0.05 are pos selection
  
  freq_table <- table(c(slac_sites, meme_sites, fel_sites))
  
  return(c(as.numeric(names(freq_table[freq_table >= 2]))))
  #return numbers which appear at least twice
  
  
}

get_z_scores <- function(recomb_proteins_pos_sel, datamonkey_codons){
  
  return(
    
    data.frame(
      
      do.call(rbind, lapply(recomb_proteins_pos_sel, function(x){
        
        unlist(zScales(x)) }
            
            ))
      
      ))
  
  #gets the Z-score for each amino acid in the list, then makes it a vector
  #then binds those vectors to a data frame and returns it
  
}

get_pos_sel_proteins <- function(recomb_proteins, datamonkey_codons){
  
  return(lapply(recomb_proteins, function(x){  
    
    return(x[datamonkey_codons])
    #get positively selected amino acids at datamonkey positions
    #returns as list of amino acids
    
  }))
  
}

get_dapc_analysis <- function(zscores, num_tests, pipeline_path){
  
  #supertype and cluster used interchangeably here
  
  cluster_estimates <- c()
  
  cluster_estimates <- c(cluster_estimates, replicate(num_tests, length(
    find.clusters(zscores,
                  n.pca=ncol(zscores), 
                  max.n.clust=40)$size)))
  
  #estimate number of clusters based on BIC graph num_tests number of times
  
  cluster_estimates_table <- table(cluster_estimates)
  
  cluster_mode <- mean(c(as.numeric(names(cluster_estimates_table[cluster_estimates_table == max(cluster_estimates_table)]))))
  #gets mode of cluster estimates; automatically rounds down if the mode is a decimal
  
  cat("\nMode Cluster Estimate: ", cluster_mode, "\n")
  cat("Median Cluster Estimate: ", median(cluster_estimates), "\n")
  cat("Average Cluster Estimate: ", mean(cluster_estimates), "\n")
  cat("Range Cluster Estimate: ", diff(range(cluster_estimates)), "\n")
  cat("Min Cluster Estimate: ", min(range(cluster_estimates)), "\n")
  cat("Max Cluster Estimate: ", max(range(cluster_estimates)), "\n")
  
  #report other summary stats of the cluster estimates
  
  cluster_data <- find.clusters(zscores, 
                                n.clust=cluster_mode,
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
  
  cat("\nFinding optimal number of PCs... This may take a moment...\n")
  
  pdf(paste0(pipeline_path,"A_optimization_graph_MGR.pdf"), width = 10, height = 10)

  a_spline_data <- optim.a.score(test_dapc, n.sim=1000)
  #finds the optimal number of principal components
  
  dev.off()
  
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
  
  recomb_by_supertype_graph <- (ggplot(indv_by_supertype_df,
               aes(x=Recombinant, 
                   y=Supertype, 
                   col=Supertype)) +
    geom_point(size = 4, alpha=0.7) +
      labs(title="Recombinant by Supertype") +
      theme(
        axis.text.x = element_text(size = 6, 
                                   angle = 90, 
                                   hjust = 1, 
                                   vjust=0.5, 
                                   margin = margin(t = 5)), 
                                   legend.position = "right",
        axis.text.y = element_text(size = 45),
        axis.title = element_text(size = 75, face = "bold"),
        legend.title = element_text(size = 60, face = "bold"),
        legend.text = element_text(size = 45),
        plot.title = element_text(size = 90, face = "bold", hjust=0.5)) +
      
      guides(color = guide_legend(ncol = 1, reverse = TRUE)))
  
  ggsave(paste0(pipeline_path, "Recombinant_by_Supertype_Graph.pdf"), 
         plot = recomb_by_supertype_graph, 
         width = 0.111*nrow(indv_by_supertype_df), 
         height = 1.1 * nlevels(indv_by_supertype_df$Supertype), 
         limitsize = FALSE)
  
  #saves the graph of the recombinants by supertype to pipeline folder
  #each supertype is a different color, and each recombinant has an associated dot
  #of that color at a height correlated with the supertype number
  
  cat("\n\nGenerating Multi-Dimensional Supertype Atlas PDF... Please wait...\n")
  
  da_columns <- paste0("LD", 1:final_dapc$n.da)
  
  all_pairs <- combn(da_columns, 2) 
  
  pdf(paste0(pipeline_path,"Recombinant_by_Supertype_Atlas.pdf"), width = 12, height = 10)
  
  ind_coords <- data.frame(final_dapc$ind.coord, 
                           Cluster = as.factor(final_dapc$assign))
  
  center_coords <- data.frame(final_dapc$grp.coord, 
                              Cluster = rownames(final_dapc$grp.coord))

  pair_list <- as.list(as.data.frame(all_pairs))
  
  invisible(lapply(pair_list, function(pair) {
    
    x_axis <- pair[1]
    y_axis <- pair[2]
    
    p <- ggplot() +
      geom_point(data = ind_coords, 
                 aes(x = .data[[x_axis]], y = .data[[y_axis]], color = Cluster), 
                 size = 2, alpha = 0.3) +
      
      geom_point(data = center_coords, 
                 aes(x = .data[[x_axis]], y = .data[[y_axis]], fill = Cluster), 
                 size = 10, shape = 21, color = "black", stroke = 1.5) +
      
      geom_text(data = center_coords, 
                aes(x = .data[[x_axis]], y = .data[[y_axis]], label = Cluster), 
                color = "white", fontface = "bold", size = 5) +
      
      theme_minimal() +
      labs(title = paste("Evolutionary Map:", x_axis, "vs", y_axis),
           x = x_axis,
           y = y_axis) +
      
      theme(
        plot.title = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        legend.position = "none" 
      )
    
    print(p) 
  }))
  
  dev.off() 
  
  cat("\nAtlas Complete!\n")
  
  #this gets all of the linear discriminants and plots them against each other in a multi-page plot 
  #helps visualize X-dimensional space, where X is the number of LDs
  
  return(final_dapc)
  
}

get_supertypes <- function(final_dapc){
  
  seahorse_names <- sub("_.*", "", names(final_dapc$grp))
  #remove underscore and all that follows to get base population name per recombinant
  
  return(data.frame(INDIVIDUAL=seahorse_names,
    SUPERTYPE = final_dapc$grp))
  
  #get population for each recombinant and its supertype (grp) 
  
}

get_glm_analyses <- function(sm_data, kept_microbe_names, used_supertypes) {
  
  cat("\nCalculating GLM... This may take a moment...\n")
  
  supertype_microbe_combos <- expand.grid(Microbe = kept_microbe_names,
                                          Supertype = used_supertypes,
                                          stringsAsFactors = FALSE)
  
  #all SupertypeMicrobe combinations for GLM analyses
  #helps establish statistical patterns between any combo
  
  GLMresults <- apply(supertype_microbe_combos, MARGIN = 1, FUN = function(supertype_microbe) {
    
    supertype <- unname(supertype_microbe["Supertype"])
    microbe   <- unname(supertype_microbe["Microbe"])
    
    glm_analysis <- suppressWarnings(coef(summary(glm(sm_data[, supertype] ~ sm_data[, microbe], 
               family = quasibinomial(link = "logit")))))
    
    #the actual GLM analysis for each association
    
    data.frame(
      SUPERTYPE_MICROBE = paste0(supertype, microbe), 
      MICROBE = microbe, 
      SUPERTYPE = supertype, 
      SLOPE = glm_analysis[2], 
      SE = glm_analysis[4], 
      WALDZ = (glm_analysis[2] / glm_analysis[4]),
      p=glm_analysis[8],
      stringsAsFactors = FALSE
    )
  })
  
  #rid of pointless rownames by unnaming them
  #run GLM analysis between all supertypes and microbes 
  #(accessing parts of expand.grid is faster than nested for loops)
  
  final_glm_df <- do.call(rbind, GLMresults)
  #all GLM analyses info bound into dataframe
  
  final_glm_df$SEQ_BONFERRONI_p <- p.adjust(final_glm_df$p, method = "holm")
  #Sequential Bonferroni correction (Holm-Bonferroni) along the entire dataframe
  
  final_glm_df$SIGNIFICANCE <- as.character(symnum(final_glm_df$SEQ_BONFERRONI_p, 
                                      corr = FALSE, 
                                      na = FALSE, 
                                      cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "")))
  
  #stars based on the significance of the Bonferroni-corrected association
  #doesnt try to reformat associations, just does 1:1 significance:symbol mapping
  #if there's an NA, prints "" insteas of R's default "?"
  
  return(final_glm_df)
  
}                       

get_waldz_heatmap <- function(pipeline_path, kept_microbe_data, glm_df, num_microbes, heatmap_title, heatmap_file_name){
  
  heatmap_microbes <- names(sort(colSums(kept_microbe_data, na.rm = TRUE), decreasing = TRUE))[1:num_microbes]
  
  heatmap_df <- glm_df[glm_df$MICROBE %in% heatmap_microbes, ]
  #of the microbes which appeared at least 100 times, these are the top 32 most prevalent ones
  #gets most relevant data and also makes heatmap more readable
  
  min_z <- min(heatmap_df$WALDZ, na.rm = TRUE)
  max_z <- max(heatmap_df$WALDZ, na.rm = TRUE)
  mid_z <- (min_z + max_z)/2
  
  sm_heatmap <- ggplot(heatmap_df, aes(x = SUPERTYPE, y = MICROBE, fill = WALDZ)) +
    geom_tile(color = "black") +
    geom_text(aes(label = SIGNIFICANCE, color = WALDZ > mid_z), size = 5, vjust = 0.7)  +
    scale_color_manual(values = c("TRUE" = "white", "FALSE" = "black"), guide = "none") +
    scale_fill_gradient2(
      low = "white", mid="gray", high = "black", 
      midpoint=mid_z,
      limits = c(min_z, max_z),
      breaks = c(min_z, max_z),
      labels = ceiling(c(min_z, max_z)), 
      name = "Wald's Z"
    ) +
    guides(fill = guide_colorbar(frame.color = "black", frame.linewidth = 0.2)) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 8)
    ) +
    labs(title = heatmap_title, x = "Supertype", y = "Microbe")
  
  suppressWarnings(ggsave(paste0(pipeline_path, heatmap_file_name, ".pdf"), 
                          plot = sm_heatmap, 
                          width = 10, 
                          height = 8))
  
  #creates heatmap of the Wald's Z for each microbe-supertype association
  #prints stars based on the significance of each sequentially-Bonferroni-corrected association
  
  
}


get_microbe_supertype_analysis <- function(pipeline_path, microbe_supertype_data, num_microbes, heatmap_title, heatmap_file_name, GLM_data_csv_file_name){
  
  kept_microbe_names <- grep(x = colnames(microbe_supertype_data), pattern = "^M",value=TRUE)
  #will be the same names as the relative data
  
  used_supertypes <- grep(x = colnames(microbe_supertype_data), pattern = "^S",value=TRUE)
  #supertypes (from JLA's DAPC analysis) associated with each seahorse 
  
  kept_microbe_data <- microbe_supertype_data[, !(colnames(microbe_supertype_data) %in% used_supertypes)]
  #removes supertype values from rows, so can get count of microbe values to keep
  
  used_supertypes_data <- apply(microbe_supertype_data[, used_supertypes], 1, function(x){sub("^S", "", names(x)[as.logical(x)])})
  #subsets the part of the microbe data that involves supertypes
  #gets the name of each supertype with a 1 (logical mask) by individual seahorse name (row)
  
  used_supertypes_df <- data.frame(INDIVIDUAL = rep(microbe_supertype_data$FISH, times=lengths(used_supertypes_data)),
                                  SUPERTYPE = as.numeric(unlist(used_supertypes_data)))
  
  #gets the supertype numbers for each seahorse and puts them into a dataframe
  
  kept_microbe_data <- kept_microbe_data[, c(TRUE, colSums(kept_microbe_data[, 2:ncol(kept_microbe_data)]) > 99)]
  #parts of kept_microbe_data where there are at least 99 reads for a given microbe
  #filter applies to all non-name columns
  #then, selected columns and first column of names are selected 
  
  rownames(kept_microbe_data) <- NULL
  kept_microbe_data <- as.matrix(column_to_rownames(kept_microbe_data, "FISH"))
  #removes the $FISH column and makes it the row names instead 
  #also gets rid of real rownames to prevent overriding issues
  #this is all numbers now, so it's a much lighter matrix
  
  kept_microbe_data_relative <- decostand(kept_microbe_data, "total")
  kept_microbe_data_relative <- cbind(kept_microbe_data_relative, microbe_supertype_data[, used_supertypes])

  #makes each read count relative to the total number of reads for that seahorse
  #bind back the supertype data for GLM analyses

  absolute_microbe_supertype_glm_df <- get_glm_analyses(microbe_supertype_data, colnames(kept_microbe_data), used_supertypes)
  
  #run GLM on absolute values for microbe_supertype counts
  #exact same results for the relative data
  
  write.csv(absolute_microbe_supertype_glm_df, paste0(pipeline_path, GLM_data_csv_file_name, ".csv"))
  #saves GLM results as a CSV file
  
  get_waldz_heatmap(pipeline_path, 
                    kept_microbe_data, 
                    absolute_microbe_supertype_glm_df, 
                    32, 
                    heatmap_title, 
                    heatmap_file_name)
  
  return(kept_microbe_data)
  
  #creates heatmap for each microbe-supertype association
  
}

get_my_microbe_supertype_data <- function(analyzed_supertypes_df, microbe_supertype_data, kept_microbe_data){
  
  analyzed_supertypes_df <- analyzed_supertypes_df[analyzed_supertypes_df$INDIVIDUAL %in% microbe_supertype_data$FISH, ]
  #gets the supertypes from my DAPC analysis for the same individuals used in JLA's analysis
  #can compare supertypes for the same individuals 
  
  analyzed_supertype_matrix <- +(table(analyzed_supertypes_df) > 0)
  #creates a logical matrix where supertypes present at least once for individuals get 1s and the rest get 0's
  
  analyzed_supertype_matrix <- analyzed_supertype_matrix[, colSums(analyzed_supertype_matrix) > 0, drop=FALSE]
  #if no individual has a supertype, exclude that supertype from the matrix; keep 2D
  
  colnames(analyzed_supertype_matrix) <- sprintf("S%02d", as.numeric(colnames(analyzed_supertype_matrix)))  
  #renames supertype number to Snumber instead; 1 becomes S1, etc
  #also pads the 2 digits with 0's, as to order supertype numbers correctly on the heatmap
  
  microbe_supertype_data <- data.frame(
    FISH = rownames(kept_microbe_data),
    as.data.frame(kept_microbe_data),
    as.data.frame(analyzed_supertype_matrix),
    stringsAsFactors = FALSE
  )
  #binds the supertype matrix and microbe counts together, just with my supertype boolean rows instead of JLA's
  #also keeps the $FISH column name, since JLA had it and since it's integrated into the pipeline
  
  return(microbe_supertype_data)
  
}

main_alleles_to_supertypes <- function(ABphasepath, pipeline_path){
  
  #with good alleles, make recombinants (a1b1, a1b2, etc)
  #get FASTA DNA sequences of alleles from PHASED alleles
  #Alpha Alleles and Beta Alleles, glue them together
  
  alpha_phasepath <- gsub("AB", "A", ABphasepath)
  beta_phasepath <- gsub("AB", "B", ABphasepath)
  
  recombs <- suppressWarnings(get_recombs(alpha_phasepath, beta_phasepath))
  
  if(!dir.exists(pipeline_path)) {
    dir.create(pipeline_path)
  }
  #creates Pipeline folder
  
  write.fasta(as.list(recombs), 
              names(recombs), 
              file.out=paste0(pipeline_path, "recombs.fasta"))
  #FASTA file of recombs
  
  #run RECCO analyses on recombinants
  
  recombs <- suppressWarnings(remove_recco(recombs, pipeline_path))
  
  write.fasta(as.list(recombs), 
              names(recombs), 
              file.out=paste0(pipeline_path, "recombs_postRECCO.fasta"))
  
  #send post-RECCO alleles to DataMonkey
  
  recomb_proteins <- get_translation(recombs)
  
  recombs <- recombs[names(recombs) %in% names(recomb_proteins)]
  #only include sequences which didn't produce a * (stop codon)
  
  write.fasta(as.list(recombs), 
              names(recombs), 
              file.out=paste0(pipeline_path, "recombs_postRECCO.fasta"))
  
  #overwrite the post-RECCO recombs file with recombs which do not produce a *
  #send these to DataMonkey
  
  datamonkey_codons <- get_datamonkey(paste0(pipeline_path, "slac.csv"), 
                                      paste0(pipeline_path, "meme.json"), 
                                      paste0(pipeline_path, "fel.json"))
  
  #get positions of amino acids under positive selection in sampled pops
  
  recomb_proteins_pos_sel <- get_pos_sel_proteins(recomb_proteins, datamonkey_codons)
  #get codons of amino acids which are under positive selection in sampled pops
  
  zscores <- get_z_scores(recomb_proteins_pos_sel)
  #makes Z-scores of each immune amino acid into a data frame
  
  final_dapc <- get_dapc_analysis(zscores, 250, pipeline_path)
  #runs cluster analysis, DAPC, and a-score optimization on z-score data for amino acids
  #returns a DAPC with the optimal number of PCs
  
  analyzed_supertypes_df <- get_supertypes(final_dapc)
  #gets supertypes of all recombinants and their population
  
  microbe_supertype_data <- read.csv(paste0(pipeline_path, "GLMOTUSTv2.csv"))
  #reads per microbe for individual seahorse, whose supertypes can be accessed with the previous dataframe
  
  kept_microbe_data <- get_microbe_supertype_analysis(pipeline_path, 
                                 microbe_supertype_data, 
                                 32, 
                                 "JLA Microbe–Supertype Associations for Microbes Present 99+ Times", 
                                 "JLA_Supertype_Microbe_WaldZ_Heatmap",
                                 "JLA_GLM_analysis_data")
  #calculates the GLM for each supertype-microbe association 
  #puts results into dataframes and generates a heatmap
  
  microbe_supertype_data <- get_my_microbe_supertype_data(analyzed_supertypes_df, microbe_supertype_data, kept_microbe_data)
  
  kept_microbe_data <- get_microbe_supertype_analysis(pipeline_path, 
                                 microbe_supertype_data, 
                                 32, 
                                 "MGR Microbe–Supertype Associations for Microbes Present 99+ Times", 
                                 "MGR_Supertype_Microbe_WaldZ_Heatmap",
                                 "MGR_GLM_analysis_data")
  
}

main_alleles_to_supertypes("PHASED/p=0.5/NO_CLONES/PHASE_AB_NOCLONES/", "Pipeline/")


