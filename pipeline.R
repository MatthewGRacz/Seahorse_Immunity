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
  
  if(sum(names(alpha_seqs) != names(beta_seqs)) == 0){cat("Good to go!\n")}
  #all sequence names are the same 
  
  names <- unique(gsub("[ab]$", "", names(alpha_seqs)))
  #get names of individuals whose sequences you're looking at
  
  a1b1 <- toupper(paste0(alpha_seqs[paste0(names, "a")], substr(beta_seqs[paste0(names, "a")], 3, 272)))
  #a = 1st allele, b = 2nd allele, per individual
  #the snipped betas have 2 extra bps ahead and 1 extra bp behind, so cut them out
  #trims to functional reading frames
  a1b2 <- toupper(paste0(alpha_seqs[paste0(names, "a")], substr(beta_seqs[paste0(names, "b")], 3, 272)))
  a2b1 <- toupper(paste0(alpha_seqs[paste0(names, "b")], substr(beta_seqs[paste0(names, "a")], 3, 272)))
  a2b2 <- toupper(paste0(alpha_seqs[paste0(names, "b")], substr(beta_seqs[paste0(names, "b")], 3, 272)))
  
  return(c(
    
    setNames(a1b1, paste0(names, "_a1b1")),
    setNames(a1b2, paste0(names, "_a1b2")),
    setNames(a2b1, paste0(names, "_a2b1")),
    setNames(a2b2, paste0(names, "_a2b2"))
    
  ))
  
}

remove_recco <- function(recombs, pipeline_phasepath){
  
  recco_phasepath <- paste0(pipeline_phasepath, "recco_results.csv")
  
  recco_names <- read.delim(recco_phasepath)$Sequence
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

get_dapc_analysis <- function(zscores, num_tests, pipeline_phasepath){
  
  #supertype and cluster used interchangeably here
  
  cluster_estimates <- c()
  
  cluster_estimates <- c(cluster_estimates, replicate(num_tests, length(
    find.clusters(zscores,
                  n.pca=ncol(zscores), 
                  max.n.clust=40)$size)))
  
  #estimate number of clusters based on BIC graph num_tests number of times
  
  cluster_estimates_table <- table(cluster_estimates)
  
  cluster_mode <- mean(c(as.numeric(names(cluster_estimates_table[cluster_estimates_table >= max(cluster_estimates_table)]))))
  
  cat("\n")
  cat("Mode Cluster Estimate: ", cluster_mode, "\n")
  
  #gets mode of the cluster estimates
  
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
  
  cat("\n")
  cat("Finding optimal number of PCs... This may take a moment...\n")

  a_spline_data <- optim.a.score(test_dapc, n.sim=10)
  #finds the optimal number of principal components
  
  final_dapc <- dapc(zscores, 
                     grp = cluster_data$grp, 
                     n.pca = a_spline_data$best,
                     n.da = 1000)
  
  #the groups are the actual separated groups using the best data, accounting for every minor detail
  #DAPC can now get the best analysis of the recombinants using only the optimized number of PCs
  #it looks at the top n.pca separations and keeps those, while the rest are deleted
  
  cat("\n")
  cat(final_dapc$var * 100, "% of Variance explained using ", final_dapc$n.pca, " PCs and ", length(final_dapc$eig), "/", length(final_dapc$pca.eig), " eigenvalues")
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
  
  ggsave(paste0(pipeline_phasepath, "Recombinant_by_Supertype_Graph.pdf"), 
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
  
  pdf(paste0(pipeline_phasepath,"Recombinant_by_Supertype_Atlas.pdf"), width = 12, height = 10)
  
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
  
  cat("Atlas Complete!\n")
  
  #this gets all of the linear discriminants and plots them against each other in a multi-page plot 
  #helps visualize X-dimensional space, where X is the number of LDs
  
  return(final_dapc)
  
}

get_supertypes <- function(final_dapc){
  
  population_names <- sub("_.*", "", names(final_dapc$grp))
  #remove underscore and all that follows to get base population name per recombinant
  
  return(data.frame(POPULATION=population_names,
    SUPERTYPE = final_dapc$grp))
  
  #get population for each recombinant and its suypertype (grp) 
  
}

main <- function(ABphasepath, pipeline_phasepath){
  
  #with good alleles, make recombinants (a1b1, a1b2, etc)
  #get FASTA DNA sequences of alleles from PHASED alleles
  #Alpha Alleles and Beta Alleles, glue them together
  
  alpha_phasepath <- gsub("AB", "A", ABphasepath)
  beta_phasepath <- gsub("AB", "B", ABphasepath)
  
  recombs <- suppressWarnings(get_recombs(alpha_phasepath, beta_phasepath))
  
  if(!dir.exists(pipeline_phasepath)) {
    dir.create(pipeline_phasepath)
  }
  #creates Pipeline folder
  
  write.fasta(as.list(recombs), 
              names(recombs), 
              file.out=paste0(pipeline_phasepath, "recombs.fasta"))
  #FASTA file of recombs
  
  #run RECCO analyses on recombinants
  
  recombs <- suppressWarnings(remove_recco(recombs, pipeline_phasepath))
  
  write.fasta(as.list(recombs), 
              names(recombs), 
              file.out=paste0(pipeline_phasepath, "recombs_postRECCO.fasta"))
  
  #send post-RECCO alleles to DataMonkey
  
  recomb_proteins <- get_translation(recombs)
  
  recombs <- recombs[names(recombs) %in% names(recomb_proteins)]
  #only include sequences which didn't produce a * (stop codon)
  
  write.fasta(as.list(recombs), 
              names(recombs), 
              file.out=paste0(pipeline_phasepath, "recombs_postRECCO.fasta"))
  
  #overwrite the post-RECCO recombs file with recombs which do not produce a *
  #send these to DataMonkey
  
  datamonkey_codons <- get_datamonkey(paste0(pipeline_phasepath, "slac.csv"), 
                                      paste0(pipeline_phasepath, "meme.json"), 
                                      paste0(pipeline_phasepath, "fel.json"))
  
  #get positions of amino acids under positive selection in sampled pops
  
  recomb_proteins_pos_sel <- get_pos_sel_proteins(recomb_proteins, datamonkey_codons)
  #get codons of amino acids which are under positive selection in sampled pops
  
  zscores <- get_z_scores(recomb_proteins_pos_sel)
  #makes Z-scores of each immune amino acid into a data frame
  
  final_dapc <- get_dapc_analysis(zscores, 1, pipeline_phasepath)
  #runs cluster analysis, DAPC, and a-score optimization on z-score data for amino acids
  #returns a DAPC with the optimal number of PCs
  
  supertype_df <- get_supertypes(final_dapc)
  #gets supertypes of all recombinants and their population
  
  microbe_supertype_data <- read.csv(paste0(pipeline_phasepath, "GLMOTUSTv2.csv"))
  #reads per microbe for individual fish, whose supertypes can be accessed with the previous dataframe
  
  supertype_df$SUPERTYPE[supertype_df$POPULATION %in% microbe_supertype_data$FISH]
  #get supertypes of fish used in microbe association
  #helps associate supertypes with certain microbes
  
  kept_microbe_names <- colnames(microbe_supertype_data)[grep(x = colnames(microbe_supertype_data), pattern = "^M",value=TRUE)]
  #will be the same names as the relative data
  
  fish_names <- row.names(microbe_supertype_data)
  
  JLA_supertypes <- colnames(microbe_supertype_data)[grep(x = colnames(microbe_supertype_data), pattern = "^S",value=TRUE)]
  
  kept_microbe_supertype_data <- microbe_supertype_data[,!(colnames(microbe_supertype_data) %in% JLA_supertypes)]
  #removes supertype values from rows, so can get count of microbe values to keep
  
  kept_microbe_data <- kept_microbe_supertype_data[colSums(kept_microbe_supertype_data[, 2:ncol(kept_microbe_supertype_data)]) > 99]
  #parts of microbe_supertype_data where there are at least 99 reads for a given microbe
  
  kept_microbe_data <- as.matrix(column_to_rownames(kept_microbe_data, "FISH"))
  #removes the $FISH column and makes it the row names instead 
  #this is all numbers now, so it's a much lighter matrix
  
  View(kept_microbe_data)
  
  kept_microbe_data_relative <- decostand(kept_microbe_data, "total")
  #makes each read count relative to the total number of reads for that seahorse
  
  View(kept_microbe_data_relative)
  
  View(supertype_df)
  
  JLA_supertypes_matrix <- matrix(row.names = fish_names,
         kept_microbe_data[JLA_supertypes])
  
  View(JLA_supertypes_matrix)
  
  
  absolute_microbe_supertype_glm <- glm()
  
  
  
  
}

main("PHASED/p=0.5/NO_CLONES/PHASE_AB_NOCLONES/", "Pipeline/")


