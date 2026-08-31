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
  recco_names <- recco_results$Sequence[recco_results$`Aln.pv` <= 0.05]
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

get_jumpers <- function(pipeline_path, zscores, num_runs, num_supertypes){

  num_runs <- num_runs + 1
  
  cat("\nRunning DAPC analyses... This may take several hours...")
  
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
  
  comparison_cols <- colnames(merged_dapc)[!colnames(merged_dapc) %in% c("RECOMBINANT", "RUN_1")]
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
  
  
  jumper_df$NOISE2 <- apply(jumper_df[, run_cols], 1, function(recomb_row){
    
    freq_table <- sort(table(as.numeric(recomb_row))/length(recomb_row), decreasing = TRUE)
    
    noise <- sum(freq_table - min(freq_table)) * (1 - max(freq_table))
    
    return(noise)
    
  })
  
  jumper_df$NOISE2 <- jumper_df$NOISE2*100/max(jumper_df$NOISE2)
  
  jumper_df$NOISE3 <- (jumper_df$NOISE * jumper_df$NOISE2)
  jumper_df$NOISE3 <- jumper_df$NOISE3*100/max(jumper_df$NOISE3)
  
  jumper_df$NOISE4 <- jumper_df$UNIQUE_STS / (jumper_df$EFFECTIVE_STS_HILL_SIMPSON * jumper_df$STABILITY)
  
  jumper_df$NOISE4 <- jumper_df$NOISE4*100/max(jumper_df$NOISE4)
  
  jumper_df$NOISE5 <- apply(jumper_df[, run_cols], 1, function(recomb_row){
    
    freq_table <- table(as.numeric(recomb_row)) / length(recomb_row)
    
    freq_table <- pmax(freq_table - (1/num_supertypes), 0)
    
    return(sum(freq_table > 0))   #num unique above chance (1/num_supertypes)
    
  })
  
  jumper_df$NOISE5 <- jumper_df$NOISE5/jumper_df$EFFECTIVE_STS_HILL_SIMPSON
  
  jumper_df$NOISE6 <- apply(jumper_df[, run_cols], 1, function(recomb_row){
    
    freq_table <- table(as.numeric(recomb_row)) / length(recomb_row)
    
    freq_table <- pmax(freq_table - (1/num_supertypes), 0)
    #all assignment proportions below random chance are removed and floored to 0
    
    freq_table <- freq_table * (1/sum(freq_table))
    #normalizes it to 100%
    
    return(sum(freq_table > 0) / (1/(sum((freq_table)^2))))
    
  })
  
  jumper_df$NOISE7 <- jumper_df$NOISE6/jumper_df$NOISE5
  
  jumper_df$NOISE5 <- jumper_df$NOISE5*100/max(jumper_df$NOISE5)
  jumper_df$NOISE6 <- jumper_df$NOISE6*100/max(jumper_df$NOISE6)
  
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
  
  write_csv(jumper_df, paste0(pipeline_path, "jumpers_dataframe.csv"))
  
  write_csv(recomb_stabilities, paste0(pipeline_path, "recombinant_stabilities.csv"))
  
  write_csv(indv_stabilities, paste0(pipeline_path, "individual_stabilities.csv"))
  
  write_csv(ari_df, paste0(pipeline_path, "ari_df.csv"))

}

#runs the 1200 DAPCs, or however many, using the Hungarian Algorithm to keep the supertype numbers consistent throughout
#also calculates stability, noise, number of unique supertypes, Simpson-Hill Index, etc per recombinant

filter_jumpers <- function(jumpers, max_noise){
  
  return(jumpers[jumpers$NOISE6 <= max_noise, ])
  
}

get_selection_site_seq <- function(recombs, codon_positions) {
  sapply(recombs, function(seq) {
    chars <- strsplit(seq, "")[[1]]
    codons <- sapply(codon_positions, function(n) get_codon(chars, n))
    paste(codons, collapse = "")
  })
}

get_indv_unique_recombs <- function(fish, recombs){
  
  return(recombs[grepl(paste0("^", fish, "_"), names(recombs))])
  #returns unique alleles for that individual
  #the reason names are ordered is so that they can be fed into FABOX, which 
  #numbers unique alleles in the order that they appear
  #JLA used this in her analysis, and the order written in her thesis (pg. 69)
  #so this is beneficial to confirm her results via running the same FABOX analysis
  #although JLA's wasn't run on the GS, but a previous version
  
}

get_captive_supertypes <- function(pipeline_path, CP_pipeline_path, recombs, CP_recombs, datamonkey_codons, zscores){

  recomb_selection_seqs <- get_selection_site_seq(recombs, datamonkey_codons)
  
  CP_recomb_selection_seqs <- get_selection_site_seq(CP_recombs, datamonkey_codons)
  
  ST_AB_CP <- read.csv(paste0(CP_pipeline_path, "recombinant_stabilities.csv"))
  ST_AB_CP <- ST_AB_CP[,c("RECOMBINANT", "INDIVIDUAL", "BEST")]
  ST_AB_CP$PROTEIN <- get_translation(CP_recomb_selection_seqs)[match(ST_AB_CP$RECOMBINANT, names(CP_recomb_selection_seqs))]
  ST_AB_CP$PROTEIN <- gsub("[^A-Z]", "", ST_AB_CP$PROTEIN)
  
  jumpers <- read.csv(paste0(pipeline_path, "jumpers_dataframe.csv"))
  
  parent_supertypes <- as.factor(jumpers$BEST[match(rownames(zscores), jumpers$RECOMBINANT)])
  
  parent_DAPC <- dapc(zscores, 
                         grp = parent_supertypes, 
                         n.pca = ncol(zscores), 
                         n.da = length(levels(parent_supertypes))-1)
  
  
  CP_prediction <- predict(parent_DAPC, 
                              newdata = get_z_scores(strsplit(ST_AB_CP$PROTEIN, "")))
  
  ST_AB_CP$TRANSLATED_BEST <- sprintf("S%02d", as.numeric(as.character(CP_prediction$assign)) - 1)
  
  return(ST_AB_CP)

}

analyze_fabox <- function(pipeline_path, jumpers){

  recomb_freqs <- read.csv(paste0(pipeline_path, "fabox_results.csv"))[, c(2, 3, 5)]
  #AB groups
  
  colnames(recomb_freqs) <- c("UNIQUE_RECOMBINANT", "FREQ", "RECOMBINANT")
  
  recomb_freqs$UNIQUE_RECOMBINANT <- sprintf("AB%03d", as.numeric(recomb_freqs$UNIQUE_RECOMBINANT))
  
  recomb_freqs <- separate_rows(recomb_freqs, RECOMBINANT, sep = "\n")
  
  jumpers_dapc <- jumpers[,c("RECOMBINANT", "BEST", "STABILITY")]
  
  jumpers_dapc <- jumpers_dapc[jumpers_dapc$RECOMBINANT %in% recomb_freqs$RECOMBINANT, ]
  
  recomb_freqs$INDIVIDUAL <- gsub("_.*", "", recomb_freqs$RECOMBINANT)
  
  recomb_freqs <- merge(recomb_freqs, jumpers_dapc, "RECOMBINANT")
  
  recomb_freqs <- cbind(recomb_freqs, SUPERTYPE = data.frame(jumpers_dapc$BEST)[jumpers_dapc$RECOMBINANT %in% recomb_freqs$RECOMBINANT, ])
  
  recomb_freqs <- cbind(recomb_freqs, data.frame(INDIVIDUAL = gsub("_.*", "", recomb_freqs$RECOMBINANT)))
  
  rownames(recomb_freqs) <- NULL
  
  return(recomb_freqs)

}

get_MST_data <- function(pipeline_path, jumpers, ordered_names){

  microbe_attribute_data <- read.csv(paste0(pipeline_path, "GLMOTUSTv2.csv"))
  microbe_attribute_data <- column_to_rownames(microbe_attribute_data, "FISH")

  JLA_matrix <- microbe_attribute_data[!grepl("^M", colnames(microbe_attribute_data))]
  JLA_matrix <- as.matrix(JLA_matrix)
  
  gut_supertypes <- jumpers[jumpers$INDIVIDUAL %in% ordered_names, ]
  #filter based on noise here!
  ST_matrix <- table(gut_supertypes$INDIVIDUAL, factor(gut_supertypes$BEST, levels = 1:17))
  colnames(ST_matrix) <- sprintf("S%02d", 0:16)
  ST_matrix <- ST_matrix[ordered_names,]
  ST_matrix <- +(ST_matrix > 0)
  ST_matrix <- unclass(ST_matrix)
  class(ST_matrix) <- "matrix"
  
  
  return(list(JLA_matrix, ST_matrix))

}

get_JLA_codons <- function(){
  
  JLA_SLAC <- read.csv('Matt_Geneious/Amena - MHC Supertype Analysis-Downloaded010721/Supertyping Final/AB Positive Selection Tests/sLACAB.csv')
  JLA_MEME <- read.csv('Matt_Geneious/Amena - MHC Supertype Analysis-Downloaded010721/Supertyping Final/AB Positive Selection Tests/MEMEAB.csv')
  JLA_FEL <- read.csv('Matt_Geneious/Amena - MHC Supertype Analysis-Downloaded010721/Supertyping Final/AB Positive Selection Tests/FELAB.csv')
  JLA_codons <- table(c(JLA_SLAC$Site[JLA_SLAC$P..dN.dS...1. <= 0.05], JLA_MEME$Site[JLA_MEME$p.value <= 0.05], JLA_FEL$Site[JLA_FEL$p.value <= 0.05]))>=2
  JLA_codons <- as.numeric(as.character(names(JLA_codons)[JLA_codons]))
  #from 2018!
  
  return(JLA_codons)
  
}




get_glm_and_heatmap_STM <- function(pipeline_path, 
                                microbe_data, 
                                JLA_microbes,
                                attribute_name, 
                                attribute_data, 
                                heatmap_title, 
                                x_scale,
                                y_scale){
  

  microbe_data <- microbe_data[, colnames(microbe_data) %in% JLA_microbes]
  #microbe data is the same for the G and H (Microbe vs Supertype and Microbe vs Unique_Recombinant)
  
  attribute_data <- attribute_data[rownames(microbe_data), , drop = FALSE]
  #puts microbe names in the same order as JLA's heatmap
  
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
  
  heatmap <- suppressWarnings(ggplot(GLMresults, aes(x = ATTRIBUTE, y = MICROBE, fill = WALDZ)) +
                                geom_tile(color = "black", size=0.4) +
                                geom_text(aes(label = SIGNIFICANCE), color = "white", size = 4, vjust = 0.7) +
                                scale_fill_gradient2(
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
  

  print(heatmap)
  
  return(GLMresults)
  
}


get_MAB_data <- function(pipeline_path, recomb_freqs){
  
  microbe_attribute_data_MAB <- read.csv(paste0(pipeline_path, "GLMOTUMH-021521.csv"))
  microbe_attribute_data_MAB <- column_to_rownames(microbe_attribute_data_MAB, "FISH")
  #microbe data is the same for the G and H (Microbe vs Supertype and Microbe vs Unique_Recombinant)
  
  JLA_attribute_data <- microbe_attribute_data_MAB[, grep(x = colnames(microbe_attribute_data_MAB), pattern = "^AB", value=TRUE)]
  JLA_attribute_freqs <- colSums(JLA_attribute_data)
  
  attribute_data <- +(table(recomb_freqs[, c("INDIVIDUAL", "UNIQUE_RECOMB")]) > 0) 
  #logical vectors of if a unique_recomb is present or not per individual
  #also, recomb_freqs already named the unique_recombs with the AB03d format, so no need for renaming
  
  attribute_freqs <- colSums(attribute_data)
  #how many individuals, not recombinants, have that unique_recomb; exactly what JLA did
  
  return(list(JLA_attribute_data, JLA_attribute_freqs, attribute_data, attribute_freqs))
  
}


get_glm_and_heatmap_MAB <- function(pipeline_path, 
                                    microbe_data, 
                                JLA_microbes,
                                attribute_name, 
                                attribute_data, 
                                attribute_freqs, 
                                heatmap_title, 
                                x_scale,
                                y_scale){
  
  microbe_data <- microbe_data[, colnames(microbe_data) %in% JLA_microbes]
  #microbe data is the same for the G and H (Microbe vs Supertype and Microbe vs Unique_Recombinant)
  
  attribute_data <- attribute_data[rownames(microbe_data), , drop = FALSE]
  #puts microbe names in the same order as JLA's heatmaps
  
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
  
  print(min_z)
  print(max_z)
  
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
  
  print(heatmap)
  
  return(GLMresults)
  
}




#main_alleles_to_supertypes("PHASED/p=0.5/CLONES/6-24-2026/AB_CLONES/", "Pipeline/AB_62426/")
#main_alleles_to_supertypes("PHASED/p=0.5/CLONES/7-27-2026/AB_CLONES/", "Pipeline/AB_72726/")
#main_alleles_to_supertypes("PHASED/p=0.5/NO_CLONES/PHASE_AB_NOCLONES/", "Pipeline/AB_NOCLONES/")


alleles_to_zscores <- function(ABphasepath, pipeline_path, is_CP){
  
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
  
  invisible(readline(prompt = "\nJust created recombinants! Upload the 'recombs.fasta' file to RECCO with 10,000 permutations, savings >=5, and\nalignment p-val<=0.05, and save it as a CSV named 'recco_results.csv' to the Pipeline folder.\nOnce done, press [Enter]!"))
  
  recombs <- suppressWarnings(remove_recco(recombs, pipeline_path))
  #run RECCO analyses on recombinants
  
  recomb_proteins <- get_translation(recombs)
  
  recombs <- recombs[names(recombs) %in% names(recomb_proteins)]
  #only include sequences which didn't produce a * (stop codon)
  
  
  write.fasta(as.list(recombs), 
              names(recombs), 
              file.out=paste0(pipeline_path, "recombs_postRECCO.fasta"))
  
  #send post-RECCO alleles to DataMonkey
  
  invisible(readline(prompt = "\nFiltered from RECCO! Send 'recombs_postRECCO.fasta.' to Datamonkey's SLAC, FEL, and MEME sites,\nand save the results in the pipeline folder as 'slac.csv', 'meme.json', and 'fel.json'.\nOnce done, press [Enter]!"))
  
  if(!is_CP){
    
    datamonkey_codons <- get_datamonkey(paste0(pipeline_path, "slac.csv"), 
                                        paste0(pipeline_path, "meme.json"), 
                                        paste0(pipeline_path, "fel.json"))
    
  }
  
  else{
    
    datamonkey_codons <- get_datamonkey(paste0(gsub("CP/", "", pipeline_path), "slac.csv"), 
                                        paste0(gsub("CP/", "", pipeline_path), "meme.json"), 
                                        paste0(gsub("CP/", "", pipeline_path), "fel.json"))
    
  }
  
  
  #get positions of amino acids under positive selection in sampled pops
  
  recomb_proteins_pos_sel <- get_pos_sel_proteins(recomb_proteins, datamonkey_codons)
  #get codons of amino acids which are under positive selection in sampled pops
  
  recomb_proteins_pos_sel <<- recomb_proteins_pos_sel
  
  zscores <- get_z_scores(recomb_proteins_pos_sel)
  View(zscores)
  
  write.csv(zscores, paste0(pipeline_path, "zscores.csv"))
  
  get_jumpers(pipeline_path, zscores, 1200, 17)
  #about 17-20 seconds per DAPC for 650ish sequences (GS dataset)
  #180-212 runs/hour, so 1200 runs is about 6 hours
  #500 runs is a bit over 2 hours
  
  return(list(zscores, recombs, datamonkey_codons))
  
}

get_shared_captives <- function(ST_AB_CP, ST_matrix){
  
  #see how many captives in $TRANSLATED_BEST share significant supertypes 
  
  return(table(factor(ST_AB_CP$TRANSLATED_BEST, levels = colnames(ST_matrix))))
  
}


#clean up heatmap functions

main <- function(pipeline_path, AB_phasepath){
  
  pipeline_path <- "Pipeline2/AB_NOCLONES/"
  AB_phasepath <- "PHASED/p=0.5/NO_CLONES/AB_NOCLONES/"
  
  AB_data <- alleles_to_zscores(AB_phasepath, pipeline_path, FALSE)
  zscores <- AB_data[[1]]
  recombs <- AB_data[[2]]
  datamonkey_codons <- AB_data[[3]]
  
  CP_pipeline_path <- paste0(pipeline_path, "CP/")
  CP_AB_phasepath <- gsub("NO_CLONES", "CP", AB_phasepath)
  
  CP_data <- alleles_to_zscores(CP_AB_phasepath, CP_pipeline_path, TRUE)
  CP_zscores <- CP_data[[1]]
  CP_recombs <- CP_data[[2]]
  
  ST_AB_CP <- get_captive_supertypes(pipeline_path, CP_pipeline_path, recombs, CP_recombs, datamonkey_codons, zscores)
  
  ordered_names <- gsub("^([A-Za-z]+)(\\d+).*$", "\\1W\\2", unlist(strsplit("SI01-15 NP01-17 NP02-17 NP03-17 NP04-17 NP05-17 CC16-17 WE01-15 WE02-15 WE03-15 WE06-17 WE07-17 WE08-17 WE09-17 WE10-17 WE11-17 AK01-15 AK12-17 AK13-17 AK14-17 AK15-17", " ")))
  #from the thesis, the names of the individuals used, in order
  
  JLA_microbes <- c('M0008', 'M0003', 'M0006', 'M0016', 'M0044', 'M0018', 'M0061', 'M0081', 'M0031', 'M0014', 'M0012', 'M0007', 'M0015', 'M0019', 'M0011', 'M0005', 'M0010', 'M0030', 'M0041', 'M0085', 'M0029', 'M0020', 'M0048', 'M0013', 'M0002', 'M0004', 'M0023', 'M0032', 'M0064', 'M0001', 'M0017', 'M0009')
  #microbes JLA used for Ch 3 heatmap, in the order she has them from bottom up (how the heatmap automatically names the y-axis)
  
  indv_unique_recombs <- lapply(ordered_names, 
                                get_indv_unique_recombs, 
                                recombs=read.fasta(paste0(pipeline_path, "recombs.fasta"), as.string=TRUE))
  #done in chapter 3, before RECCO filtering and after TA filtering, hence why recombs and not recombs_postRECCO is used
  
  write.fasta(as.list(unlist(indv_unique_recombs)), 
              names=names(unlist(indv_unique_recombs)), 
              paste0(pipeline_path, "forFABOX.fasta"))
  
  #FaBox exports results as HTML, convert to CSV and use that file here
  
  jumpers <- read.csv(paste0(pipeline_path, "jumpers_dataframe.csv"))
  
  microbe_attribute_data <- read.csv(paste0(pipeline_path, "GLMOTUSTv2.csv"))
  microbe_attribute_data <- column_to_rownames(microbe_attribute_data, "FISH")
  microbe_data <- microbe_attribute_data[!grepl("^S", colnames(microbe_attribute_data))]
  
  invisible(readline(prompt = "\nAlmost done! Send the file 'forFABOX.fasta' to FaBox and export its results as an HTML file.\nSend that file to a website to convert it to CSV, then save that file in the pipeline folder\nas 'fabox_results.csv'. Once done, press [Enter]!"))
  
  recomb_freqs <- analyze_fabox(pipeline_path, jumpers)
  
  JLA_codons <- get_JLA_codons()
  #use instead of datamonkey_codons for any analysis, preferably on the GS, for comparing results
  
  jumpers <- filter_jumpers(jumpers, 95)
  #filter based on noise, stability, etc
  #max_noise set to 95
  
  MST_prep_data <- get_MST_data(pipeline_path, jumpers, ordered_names)
  
  JLA_matrix <- MST_prep_data[[1]]
  ST_matrix <- MST_prep_data[[2]]
  
  MST_data <- get_glm_and_heatmap_STM(pipeline_path, 
                                      microbe_data, 
                                      JLA_microbes,
                                      "Supertype", 
                                      ST_matrix,
                                      "Microbe–Supertype Associations", 
                                      0.3,
                                      0.3)
  
  
  JLA_MST_data <- get_glm_and_heatmap_STM(pipeline_path, 
                                          microbe_data, 
                                          JLA_microbes,
                                          "Supertype", 
                                          JLA_matrix,
                                          "JLA Microbe–Supertype Associations", 
                                          0.3,
                                          0.3)
  
  MAB_prep_data <- get_MAB_data(pipeline_path, recomb_freqs)
  
  JLA_attribute_data <- MAB_prep_data[[1]]
  JLA_attribute_freqs <- MAB_prep_data[[2]]
  attribute_data <- MAB_prep_data[[3]]
  attribute_freqs <- MAB_prep_data[[4]]
  
  
  MAB_data <- get_glm_and_heatmap_MAB(pipeline_path, 
                                      microbe_data, 
                                  JLA_microbes, 
                                  "Unique Recombinant", 
                                  attribute_data, 
                                  attribute_freqs, 
                                  "Microbe–AB Group Associations for JLA's 32 Microbes", 
                                  0.3,
                                  0.3)
  
  JLA_MAB_data <- get_glm_and_heatmap_MAB(pipeline_path, 
                                          microbe_data, 
                                      JLA_microbes, 
                                      "Unique Recombinant", 
                                      JLA_attribute_data, 
                                      JLA_attribute_freqs,
                                      "JLA Microbe–AB Group Associations for JLA's 32 Microbes", 
                                      0.3,
                                      0.3)
  
  captives_per_supertype <- get_shared_captives(ST_AB_CP, ST_matrix)
  
  
}

main("Pipeline2/AB_NOCLONES/", "PHASED/p=0.5/NO_CLONES/AB_NOCLONES/")

