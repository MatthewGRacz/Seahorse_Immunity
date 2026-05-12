library("dplyr")
library("vegan")
library("readxl")
library("adegenet")
library("BiocManager")
library("devtools")
library("rBLAST")
library("biomartr")
library("Biostrings")
library("GenomicRanges")
library("IRanges")
library("S4Vectors")
library("seqinr")
library("utils")
library("metagenomeFeatures")
library("qiime2R")
library("rentrez")
library("DBI")
library("XML")
library("adegenet")
#remotes::install_version("Peptides", version = 2.4) #version used in analysis
library("Peptides")
library("modeest")
library("umap")
library("ggplot2")
library("ggVennDiagram")
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
  
  return(data.frame(
    do.call(rbind, 
            zScales(recomb_proteins_pos_sel)), 
    row.names = names(recomb_proteins_pos_sel)))
  
}

get_pos_sel_proteins <- function(recomb_proteins, datamonkey_codons){
  
  return(c(lapply(recomb_proteins, function(x){
    
    return(paste(unlist(x)[datamonkey_codons], collapse=""))
    #get positively selected amino acids at datamonkey positions, then return as strings
    
  })))
  
}

get_dapc_analysis <- function(zscores){
  
  cluster_data <- find.clusters(zscores, max.n.clust=40, n.da=1000, n.pca=5)
  #shows BIC graph of clusters, to find optimal number of clusters
  #keep PCA=5 for there are 5 Z-scores total, even if 3 carry 90% of the variance
  #n.da set high will automatically find the maximum number of axes for the cluster comparisons
  
  dapc_data <- dapc(zscores, cluster_data$grp, n.pca=5)
  #DAPC analysis on Z-scores with the number of groups from the cluster analysis
  #still all 5 PCs kept from before for all 5 Z-scores
  
  optim.a.score(dapc_data)
  #find the optimal number of principal components
  
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
  #makes Z-scores of each immune protein into a data frame
  
  
  
}

main("PHASED/p=0.5/NO_CLONES/PHASE_AB_NOCLONES/", "Pipeline/")


