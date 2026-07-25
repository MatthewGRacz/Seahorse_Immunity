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



get_A_B_recombs <- function(alpha_phasepath, beta_phasepath){
  
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

get_AB_recombs <- function(AB_phasepath){
  
  AB_phasepath <- paste0(AB_phasepath, "phased.fasta")
  
  AB_alleles <- read.fasta(AB_phasepath, as.string = TRUE)
  
  AB_alleles <- AB_alleles[!duplicated(names(AB_alleles), fromLast = TRUE)]
  #KEY! This only looks for the newest AB alleles for each individual, so the cloned known alleles 
  #are used, and the older, original, potentially-wrong alleles are not used
  
  indv_names <- unique(gsub("[ab]$", "", names(AB_alleles)))
  #get names of individuals whose sequences you're looking at
  
  a1b1 <- toupper(paste0(substr(AB_alleles[paste0(indv_names, "a")], 1, 246), substr(AB_alleles[paste0(indv_names, "a")], 249, 518)))
  #a = 1st allele, b = 2nd allele, per individual
  #the snipped betas have 2 extra bps ahead and 1 extra bp behind, so cut them out
  #trims to functional reading frames
  a1b2 <- toupper(paste0(substr(AB_alleles[paste0(indv_names, "a")], 1, 246), substr(AB_alleles[paste0(indv_names, "b")], 249, 518)))
  a2b1 <- toupper(paste0(substr(AB_alleles[paste0(indv_names, "b")], 1, 246), substr(AB_alleles[paste0(indv_names, "a")], 249, 518)))
  a2b2 <- toupper(paste0(substr(AB_alleles[paste0(indv_names, "b")], 1, 246), substr(AB_alleles[paste0(indv_names, "b")], 249, 518)))
  
  return(c(
    
    setNames(a1b1, paste0(indv_names, "_a1b1")),
    setNames(a1b2, paste0(indv_names, "_a1b2")),
    setNames(a2b1, paste0(indv_names, "_a2b1")),
    setNames(a2b2, paste0(indv_names, "_a2b2"))
    
  ))
  
  
}

remove_recco <- function(recombs, pipeline_path){
  
  recco_path <- paste0(pipeline_path, "recco_results.csv")
  
  recco_results <- read.delim(recco_path)
  recco_names <- recco_results$Sequence[recco_results$Aln.pv <= 0.05]
  #RECCO's output file is TSV, not CSV like it says
  #where flagged recombinants with p <= 0.05 for sequence alignment are removed
  
  return(recombs[!names(recombs) %in% recco_names])
  
}

get_translation <- function(recombs){
  
  recomb_proteins <- lapply(recombs, function(x) {seqinr::translate(seqinr::s2c(x))})
  
  return(recomb_proteins[!sapply(recomb_proteins, function(x) "*" %in% x)])
  #if any have a stop codon by chance of such bps being concatenated, filter them out
  
}

get_datamonkey <- function(slac, meme, fel){
  
  slac <- read.csv(slac)
  
  slac_sites <- c(slac[slac$"P..dN.dS...1." <= 0.05, ]$"Site")
  #CSV compresses headers; this corresponds with the positive selection header
  #P [dN/dS > 1]
  
  meme <- fromJSON(meme)$MLE$content$`0`
  #where the data lies
  
  meme_sites <- c(which(meme[,7] <= 0.05))
  #7th column is p values, and all p values < 0.05 are pos selection
  
  fel <- fromJSON(fel)$MLE$content$`0`
  
  fel_sites <- c(which(fel[,5] <= 0.05))
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

  a_spline_data <- optim.a.score(test_dapc, n.sim=100)
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
    
    glm_analysis <- suppressWarnings(coef(summary(glm(attribute_data[, attribute] ~ microbe_data[, microbe], 
                                                      family = quasibinomial(link = "logit")))))
    
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
    scale_color_manual(values = c("TRUE" = "white", "FALSE" = "black"), guide = "none") +
    scale_fill_gradient2(
      low = "white", mid="#A6A6A6", high = "black", 
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
      axis.text.x = element_text(size = x_scale * 13), 
      axis.text.y = element_text(size = y_scale * 30), 
      plot.title = element_text(size = 70 * x_scale, face = "bold"), 
      axis.title.x = element_text(size = x_scale * 50, face = "bold"),
      axis.title.y = element_text(size = y_scale * 50, face = "bold")
    ) +
    labs(title = heatmap_title, x = attribute_name , y = "Microbe"))
  
  #suppressWarnings(ggsave(paste0(pipeline_path, heatmap_file_name, ".pdf"), plot = heatmap, width = x_scale * length(unique(GLMresults$ATTRIBUTE)), height = num_microbes * y_scale, limitsize = FALSE))
  
  print(heatmap)
  
  return(GLMresults)
  
}

main_alleles_to_supertypes <- function(ABphasepath, pipeline_path){
  
  ABphasepath <- "PHASED/p=0.5/CLONES/6-24-2026/AB_CLONES/"
  pipeline_path <- "Pipeline/"
  
  #with good alleles, make recombinants (a1b1, a1b2, etc)
  #get FASTA DNA sequences of alleles from PHASED alleles
  #Alpha Alleles and Beta Alleles, glue them together
  
  alpha_phasepath <- gsub("AB", "A", ABphasepath)
  beta_phasepath <- gsub("AB", "B", ABphasepath)
  
  recombs <- suppressWarnings(get_AB_recombs(ABphasepath))
  
  if(!dir.exists(pipeline_path)) {
    dir.create(pipeline_path)
  }
  #creates Pipeline folder
  
  recombs <- recombs[!grepl("^TA", names(recombs))]
  #removes TA individuals before ch 4 analyses such as RECCO and DataMonkey
  
  
  write.fasta(as.list(recombs), 
              names(recombs), 
              file.out=paste0(pipeline_path, "recombs.fasta"))
  #FASTA file of recombs
  
  invisible(readline(prompt = "\nJust created recombinants! Upload the 'recombs.fasta' file to RECCO with 10,000 permutations, savings >=5, and\nalignment p-val<=0.05, and save it as a CSV named 'recco.results.csv' to the Pipeline folder.\nOnce done, press [Enter]!"))
  
  recombs <- suppressWarnings(remove_recco(recombs, pipeline_path))
  #run RECCO analyses on recombinants
  
  recomb_proteins <- get_translation(recombs)
  
  recomb_proteins <<- recomb_proteins
  
  recombs <- recombs[names(recombs) %in% names(recomb_proteins)]
  #only include sequences which didn't produce a * (stop codon)
  
  
  write.fasta(as.list(recombs), 
              names(recombs), 
              file.out=paste0(pipeline_path, "recombs_postRECCO.fasta"))
  
  #send post-RECCO alleles to DataMonkey
  
  recombs <<- recombs
  
  invisible(readline(prompt = "\nFiltered from RECCO! Send 'recombs_postRECCO.fasta.' to Datamonkey's SLAC, FEL, and MEME sites,\nand save the results in the pipeline folder as 'slac.csv', 'meme.json', and 'fel.json'.\nOnce done, press [Enter]!"))
  
  datamonkey_codons <- get_datamonkey(paste0(pipeline_path, "slac.csv"), 
                                      paste0(pipeline_path, "meme.json"), 
                                      paste0(pipeline_path, "fel.json"))
  
  #get positions of amino acids under positive selection in sampled pops
  
  recomb_proteins_pos_sel <- get_pos_sel_proteins(recomb_proteins, datamonkey_codons)
  #get codons of amino acids which are under positive selection in sampled pops
  
  recomb_proteins_pos_sel <<- recomb_proteins_pos_sel
  
  zscores <- get_z_scores(recomb_proteins_pos_sel)
  zscores <<- zscores
  View(zscores)
  #makes Z-scores of each immune amino acid into a data frame
  
  final_dapc <- get_dapc_analysis(zscores, 1, pipeline_path)
  final_dapc <<- final_dapc
  #runs cluster analysis, DAPC, and a-score optimization on z-score data for amino acids
  #returns a DAPC with the optimal number of PCs
  
  analyzed_supertypes_df <- get_supertypes(final_dapc)
  analyzed_supertypes_df <<- analyzed_supertypes_df
  #gets supertypes of all recombinants and their population
  
  View(analyzed_supertypes_df)
  
  View(table(final_dapc$grp))

  get_indv_unique_recombs <- function(fish, recombs){

    return(recombs[grepl(paste0("^", fish, "_"), names(recombs))])
    #returns unique alleles for that individual
    #the reason names are ordered is so that they can be fed into FABOX, which 
    #numbers unique alleles in the order that they appear
    #JLA used this in her analysis, and the order written in her thesis (pg. 69)
    #so this is beneficial to confirm her results via running the same FABOX analysis
    #although JLA's wasn't run on the GS, but a previous version
      
  }
  
  microbe_attribute_data <- read.csv(paste0(pipeline_path, "GLMOTUSTv2.csv"))
  microbe_attribute_data <- column_to_rownames(microbe_attribute_data, "FISH")
  microbe_data <- microbe_attribute_data[, grep(x = colnames(microbe_attribute_data), pattern = "^M", value=TRUE)]
  
  ordered_names <- gsub("^([A-Za-z]+)(\\d+).*$", "\\1W\\2", unlist(strsplit("SI01-15 NP01-17 NP02-17 NP03-17 NP04-17 NP05-17 CC16-17 WE01-15 WE02-15 WE03-15 WE06-17 WE07-17 WE08-17 WE09-17 WE10-17 WE11-17 AK01-15 AK12-17 AK13-17 AK14-17 AK15-17", " ")))
  
  indv_unique_recombs <- lapply(ordered_names, get_indv_unique_recombs, recombs=read.fasta(paste0(pipeline_path, "recombs.fasta"), as.string=TRUE))
  #done before RECCO filtering and after TA filtering
  
  write.fasta(as.list(unlist(indv_unique_recombs)), names=names(unlist(indv_unique_recombs)), paste0(pipeline_path, "forFABOX.fasta"))
  
  #print(length(unique(unlist(indv_unique_recombs))))
  
  #FaBox exports results as HTML, convert to CSV and use that file here
  
  invisible(readline(prompt = "\nAlmost done! Send the file 'forFABOX.fasta' to FaBox and export its results as an HTML file.\nSend that file to a website to convert it to CSV, then save that file in the pipeline folder\nas 'fabox_results.csv'. Once done, press [Enter]!"))
  
  recomb_freqs <- read.csv(paste0(pipeline_path, "fabox_results.csv"))[, c(2, 3, 5)]
  
  colnames(recomb_freqs) <- c("UNIQUE_RECOMB", "FREQ", "RECOMB")
  
  recomb_freqs$UNIQUE_RECOMB <- sprintf("AB%03d", as.numeric(recomb_freqs$UNIQUE_RECOMB))
  
  recomb_freqs <- separate_rows(recomb_freqs, RECOMB, sep = "\n")
  
  recomb_freqs <- cbind(recomb_freqs, SUPERTYPE = data.frame(final_dapc$grp)[c(recomb_freqs$RECOMB) , ])
  
  recomb_freqs <- cbind(recomb_freqs, data.frame(INDIVIDUAL = gsub("_.*", "", recomb_freqs$RECOMB)))
  
  recomb_freqs <- cbind(recomb_freqs, microbe_data[recomb_freqs$INDIVIDUAL, ])
  
  rownames(recomb_freqs) <- NULL
  
  recomb_freqs <- recomb_freqs
  
  #View(recomb_freqs)
  
  JLA_microbes <- c('M0008', 'M0003', 'M0006', 'M0016', 'M0044', 'M0018', 'M0061', 'M0081', 'M0031', 'M0014', 'M0012', 'M0007', 'M0015', 'M0019', 'M0011', 'M0005', 'M0010', 'M0030', 'M0041', 'M0085', 'M0029', 'M0020', 'M0048', 'M0013', 'M0002', 'M0004', 'M0023', 'M0032', 'M0064', 'M0001', 'M0017', 'M0009')
  #microbes JLA used for Ch 3 heatmap, in the order she has them from bottom up (how the heatmap automatically names the y-axis)
  
  microbe_attribute_data <- read.csv(paste0(pipeline_path, "GLMOTUSTv2.csv"))
  microbe_attribute_data <- column_to_rownames(microbe_attribute_data, "FISH")
  microbe_data <- microbe_attribute_data[, colnames(microbe_attribute_data) %in% JLA_microbes]
  #microbe data is the same for the G and H (Microbe vs Supertype and Microbe vs Unique_Recombinant)
  
  JLA_MST_data <- get_glm_and_heatmap(pipeline_path, 
                      microbe_attribute_data, 
                      JLA_microbes,
                      "S", 
                      "Supertype", 
                      NULL, 
                      NULL, 
                      "JLA Microbe–Supertype Associations for Microbes Present 99+ Times", 
                      "JLA_Supertype_Microbe_WaldZ_Heatmap",
                      0.3,
                      0.3)
  
  attribute_data <- +(table(analyzed_supertypes_df) > 0) 
  attribute_data <- attribute_data[rownames(attribute_data) %in% rownames(microbe_data), ]
  colnames(attribute_data) <- sprintf("S%02d", as.numeric(colnames(attribute_data)))
  #logical vectors of if a supertype is present or not per individual, just with my supertypes
  #only for individuals whose microbes were analyzed
  #renames their supertype numbers to the S02d format (eg. S01, S23)
  
  attribute_freqs <- table(sprintf("S%02d", as.numeric(final_dapc$grp)))
  #the frequency of each supertype from my analysis
  #also renames final_dapc's supertype numbers to the S02d format (eg. S01, S23)
  
  #MGR Microbe vs Supertype Analysis
  
  MGR_MST_data <- get_glm_and_heatmap(pipeline_path, 
                      microbe_attribute_data, 
                      JLA_microbes, 
                      "S", 
                      "Supertype", 
                      attribute_data, 
                      attribute_freqs, 
                      "MGR Microbe–Supertype Associations for Microbes Present 99+ Times", 
                      "MGR_Supertype_Microbe_WaldZ_Heatmap",
                      0.3,
                      0.3)
  
  
  microbe_attribute_data <- read.csv(paste0(pipeline_path, "GLMOTUMH-021521.csv"))
  microbe_attribute_data <- column_to_rownames(microbe_attribute_data, "FISH")
  #microbe data is the same for the G and H (Microbe vs Supertype and Microbe vs Unique_Recombinant)
  
  attribute_data <- microbe_attribute_data[, grep(x = colnames(microbe_attribute_data), pattern = "^AB", value=TRUE)]
  attribute_freqs <- colSums(attribute_data)
  
  JLA_MAB_data <- get_glm_and_heatmap(pipeline_path, 
                      microbe_attribute_data, 
                      JLA_microbes, 
                      "AB", 
                      "Unique Recombinant", 
                      NULL, 
                      attribute_freqs,
                      "JLA Microbe–Unique Recombinant Associations for All Microbes", 
                      "JLA_UniqueRecombinant_Microbe_WaldZ_Heatmap",
                      0.3,
                      0.3)
  
  attribute_data <- +(table(recomb_freqs[, c("INDIVIDUAL", "UNIQUE_RECOMB")]) > 0) 
  #logical vectors of if a unique_recomb is present or not per individual
  #also, recomb_freqs already named the unique_recombs with the AB03d format, so no need for renaming
  
  attribute_freqs <- colSums(attribute_data)
  #how many individuals, not recombinants, have that unique_recomb; exactly what JLA did
  
  View(recomb_freqs)
  
  MGR_MAB_data <- get_glm_and_heatmap(pipeline_path, 
                      microbe_attribute_data, 
                      JLA_microbes, 
                      "AB", 
                      "Unique Recombinant", 
                      attribute_data, 
                      attribute_freqs, 
                      "MGR Microbe–Unique Recombinant Associations for All Microbes", 
                      "MGR_UniqueRecombinant_Microbe_WaldZ_Heatmap",
                      0.3,
                      0.3)
  
  #get the glm data from the microbe-supertype run 
  #get the supertype of the significant associations (where significance != "")
  #get the microbe of the significant associations
  #get the glm data from the microbe-unique_recomb run
  #get the unique_recombs of the significant associations (where significance != "")
  #get the supertypes of individuals with that unique_recomb (from recomb_freqs)
  #get the microbe of the significant associations
  
  #printing all of this data out, combined with PHASE of the captives and comparison 
  #of their supertypes with these supertypes, allows us to see if JLA had good judgement 
  #in her chapter 4 experiment of which microbes she picked and for what supertype of captive seahorse
  
  microbe_taxonomy_data <- read_xlsx(paste0(pipeline_path, "GLMOTUMH_output_MS2FIG6.xlsx"), sheet="p-value for Fig6")
  
  microbe_taxonomy_data <- na.omit(data.frame(invisible(cbind(MICROBE = microbe_taxonomy_data$OTUs, 
                                                              PHYLUM = microbe_taxonomy_data$Phylum, 
                                                              ORDER = microbe_taxonomy_data$Order, 
                                                              GENUS = microbe_taxonomy_data$Genus))))
  
  
  JLA_MST_data <- JLA_MST_data[JLA_MST_data$SIGNIFICANCE != "", ]
  JLA_MAB_data <- JLA_MAB_data[JLA_MAB_data$SIGNIFICANCE != "", ]
  MGR_MST_data <- MGR_MST_data[MGR_MST_data$SIGNIFICANCE != "", ]
  MGR_MAB_data <- MGR_MAB_data[MGR_MAB_data$SIGNIFICANCE != "", ]
  
  JLA_MST_df <- data.frame(
    MICROBE = JLA_MST_data$MICROBE,
    GENUS = microbe_taxonomy_data$GENUS[match(JLA_MST_data$MICROBE, microbe_taxonomy_data$MICROBE)],
    SUPERTYPE = as.numeric(gsub("S", "", JLA_MST_data$ATTRIBUTE))
  )
  
  JLA_MAB_df <- data.frame(
    MICROBE = JLA_MAB_data$MICROBE,
    GENUS = microbe_taxonomy_data$GENUS[match(JLA_MAB_data$MICROBE, microbe_taxonomy_data$MICROBE)],
    UNIQUE_RECOMB = JLA_MAB_data$ATTRIBUTE,
    SUPERTYPE = recomb_freqs$SUPERTYPE[match(JLA_MAB_data$ATTRIBUTE, recomb_freqs$UNIQUE_RECOMB)]
  )
  
  MGR_MST_df <- data.frame(
    MICROBE = MGR_MST_data$MICROBE,
    GENUS = microbe_taxonomy_data$GENUS[match(MGR_MST_data$MICROBE, microbe_taxonomy_data$MICROBE)],
    SUPERTYPE = as.numeric(gsub("S", "", MGR_MST_data$ATTRIBUTE))
  )
  
  MGR_MAB_df <- data.frame(
    MICROBE = MGR_MAB_data$MICROBE,
    GENUS = microbe_taxonomy_data$GENUS[match(MGR_MAB_data$MICROBE, microbe_taxonomy_data$MICROBE)],
    UNIQUE_RECOMB = MGR_MAB_data$ATTRIBUTE,
    SUPERTYPE = recomb_freqs$SUPERTYPE[match(MGR_MAB_data$ATTRIBUTE, recomb_freqs$UNIQUE_RECOMB)]
  )
  
  View(JLA_MST_df)
  View(JLA_MAB_df)
  
  supertype_mapping <- merge(MGR_MST_df, JLA_MST_df, 
                             by = "MICROBE", 
                             suffixes = c("_MGR", "_JLA"))
  
  View(supertype_mapping)
  
  
  

  View(MGR_MST_df)
  View(MGR_MAB_df)
  

  
}

main_alleles_to_supertypes("PHASED/p=0.5/CLONES/6-24-2026/AB_CLONES/", "Pipeline/")


