library("Biostrings")
library("seqinr")
library("utils")
library("adegenet")
#remotes::install_version("Peptides", version = 2.4) #version used in analysis
library("Peptides")
library("ggplot2")
library("jsonlite")

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

get_dapc_analysis <- function(zscores, num_tests){
  
  cluster_estimates <- c()
  
  cluster_estimates <- c(cluster_estimates, replicate(num_tests, length(find.clusters(zscores,
                          n.pca=ncol(zscores), max.n.clust=50)$size)))
  
  #estimate number of clusters based on BIC graph num_tests number of times
  
  cluster_estimates_table <- table(cluster_estimates)
  
  cluster_mode <- mean(c(as.numeric(names(cluster_estimates_table[cluster_estimates_table >= max(cluster_estimates_table)]))))
  
  #gets mode of the cluster estimates
  
  cat("Mode Cluster Estimate: ", cluster_mode, "\n")
  cat("Median Cluster Estimate: ", median(cluster_estimates), "\n")
  cat("Average Cluster Estimate: ", mean(cluster_estimates), "\n")
  cat("Range Cluster Estimate: ", diff(range(cluster_estimates)), "\n")
  cat("Min Cluster Estimate: ", min(range(cluster_estimates)), "\n")
  cat("Max Cluster Estimate: ", max(range(cluster_estimates)), "\n")
  
  
  cluster_data <- find.clusters(zscores, 
                                n.clust=cluster_mode,
                                n.pca=ncol(zscores))
  #shows BIC graph of clusters, to find optimal number of clusters
  #keep n.pca=ncol(zscores) because you want to use all the variance possible, to overshoot
  #the mode cluster number is used as the number of clusters
  
  test_dapc <- dapc(zscores, 
                    cluster_data$grp, 
                    n.pca=ncol(zscores), 
                    n.da=ncol(zscores))

  #all PCs kept from find.clusters for all amino acids' Z-scores
  #overshooting helps use all values between 1 and the overshot npca number
  #helps to find the most conservative, statistically probable npca number
  #n.da maxes out at n.pca-1, finds the maximum number of axes for the cluster comparisons; 
  #setting n.da too high will cause it to self-correct to the maximum possible value

  a_spline_data <- optim.a.score(test_dapc)
  #finds the optimal number of principal components
  
  final_dapc <- dapc(zscores, 
                     grp = cluster_data$grp, 
                     n.pca = a_spline_data$best,
                     n.da = a_spline_data$best)
  
  #the groups are the actual separated groups using the best data, accounting for every minor detail
  #DAPC can now get the best analysis of the recombinants using only the optimized number of PCs
  #it looks at the top n.pca separations and keeps those, while the rest are deleted
  
  return(final_dapc)
  
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
  
  best_dapc <- get_dapc_analysis(zscores, 3)
  #runs cluster analysis, DAPC, and a-score optimization on z-score data for amino acids
  #returns a DAPC with the optimal number of PCs
  
  print(best_dapc)
  
  
}

main("PHASED/p=0.5/NO_CLONES/PHASE_AB_NOCLONES/", "Pipeline/")


