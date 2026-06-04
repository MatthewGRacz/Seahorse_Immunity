library("seqinr")
library("waldo")

setwd("/Users/mattracz/Projects/Wilson_Lab")

analyze_out <- function(outpath, pipeline_path, p){
    
  outfile <- readLines(outpath)
  
  type_seq <- sub(".*PHASE_([AB]+).*", "\\1", outpath)
  #grabs A, B, or AB in phasepath name to classify allele as Alpha, Beta, or AlphaBeta
  
  #get name of allele
  #get allelic sequence
  #get frequency of alleles
  #get number of flagged loci per allele
  #get length of allele
  #get dataframe (as CSV) of confidences at each loci
  
  alleles_by_indv <- outfile[(grep("BEGIN BESTPAIRS_SUMMARY", outfile, value=FALSE)+1):(grep("END BESTPAIRS_SUMMARY", outfile, value=FALSE)-1)]
  
  alleles_by_indv <- strsplit(gsub(":", ",", gsub("\\(|\\)", "", gsub(" ", "", alleles_by_indv))), ",")
  #[1] is name, [2] is allele1 number, [3] is allele2 number
  
  alleles_by_freq <- outfile[(grep("BEGIN LIST_SUMMARY", outfile, value=FALSE)+1):(grep("END LIST_SUMMARY", outfile, value=FALSE)-1)]
  
  phased_allele_freqs <- as.integer(sub(".*\\s+([0-9]+).*$", "\\1", alleles_by_freq))
  
  #where says "positions of loci," get the loci positions
  
  allele_probs <- outfile[(grep("BEGIN PHASEPROBS", outfile, value=FALSE)+1):(grep("END PHASEPROBS", outfile, value=FALSE)-1)]
  
  allele_probs <- gsub("=", "1.00", allele_probs)
  
  allele_probs <- do.call(rbind, lapply(allele_probs, function(x){
    
    return(as.numeric(unlist(strsplit(x, "\\s+"))))
    
  }))
  
  loci_positions <- outfile[(grep("Positions of loci: ", outfile, value=FALSE))]
  loci_positions <- regmatches(loci_positions, gregexpr("[0-9]+", loci_positions))[[1]]
  
  colnames(allele_probs) <- loci_positions
  
  low_conf_loci <- apply(allele_probs, 1, function(x){
    
    loci <- names(x[x < p])
    
    if(length(loci) == 0){
      return("0") 
    }
    
    return(paste(unname(loci), collapse = ", "))
    
  })
  
  input_alleles_data <- data.frame(
    
    INDIVIDUAL = sapply(alleles_by_indv, `[`, 1),
    NUM_LOW_CONFIDENCE = rowSums(allele_probs < p),
    LOW_CONFIDENCE_LOCI = low_conf_loci,
    NUM_LOCI = length(loci_positions),
    allele_probs,
    
    stringsAsFactors = FALSE,
    check.names = FALSE
    
  )
  
  #View(input_alleles_data)
  
  write.csv(input_alleles_data, paste0(pipeline_path, "input_", type_seq, "_alleles_data.csv"))
  
  phased_allele_seqs <- suppressWarnings(toupper(read.fasta(sub("seqphase.out", "phased.fasta", outpath), as.string=TRUE)))
  #get allelic sequences
  
  phased_allele_lengths <- nchar(phased_allele_seqs)
  
  phased_individuals <- rep(sapply(alleles_by_indv, `[`, 1), each = 2)
  
  phased_allele_numbers <- c(rbind(as.integer(sapply(alleles_by_indv, `[`, 2)), as.integer(sapply(alleles_by_indv, `[`, 3))))
  
  phased_allele_data <- data.frame(
    
    INDIVIDUAL = phased_individuals,
    SEQUENCE = phased_allele_seqs[phased_allele_numbers],
    LENGTH = phased_allele_lengths[phased_allele_numbers],
    ALLELE_NUM = phased_allele_numbers,
    FREQUENCY = phased_allele_freqs[phased_allele_numbers],
    
    stringsAsFactors = FALSE
    
  )
  
  View(phased_allele_data)
  
  #write.csv(phased_allele_data, paste0(pipeline_path, "phased_", type_seq, "_alleles_data.csv"))
  
  p_phased <- outfile[(grep("BEGIN COMMAND_LINE ", outfile, value=FALSE)+1):(grep("END COMMAND_LINE", outfile, value=FALSE)-1)]
  
  p_phased <- as.numeric(gsub(".*\\s+-p([0-9]+\\.[0-9]+)\\s.*", "\\1", p_phased))
  
  return(data.frame(
    
    FILE = sub("/seqphase.out", "", outpath),
    P = p_phased,
    TYPE = type_seq,
    NUM_INDIVIDUALS = length(alleles_by_indv),
    NUM_UNIQUE_ALLELES = max(phased_allele_numbers),
    NUM_LOCI = length(loci_positions),
    MAX_LOW_CONFIDENCE_LOCI = max(input_alleles_data$NUM_LOW_CONFIDENCE)
  
  ))
  
}

get_recombs <- function(alpha_seqs, beta_seqs){
  
  indv_names <- unique(gsub("[ab]$", "", names(alpha_seqs)))
  #get names of individuals whose sequences you're looking at
  
  a1b1 <- toupper(paste0(alpha_seqs[paste0(indv_names, "a")], beta_seqs[paste0(indv_names, "a")]))
  #a = 1st allele, b = 2nd allele, per individual
  #the snipped betas have 2 extra bps ahead and 1 extra bp behind, so cut them out
  #trims to functional reading frames
  a1b2 <- toupper(paste0(alpha_seqs[paste0(indv_names, "a")], beta_seqs[paste0(indv_names, "b")]))
  a2b1 <- toupper(paste0(alpha_seqs[paste0(indv_names, "b")], beta_seqs[paste0(indv_names, "a")]))
  a2b2 <- toupper(paste0(alpha_seqs[paste0(indv_names, "b")], beta_seqs[paste0(indv_names, "b")]))
  
  return(c(
    
    setNames(a1b1, paste0(indv_names, "_a1b1")),
    setNames(a1b2, paste0(indv_names, "_a1b2")),
    setNames(a2b1, paste0(indv_names, "_a2b1")),
    setNames(a2b2, paste0(indv_names, "_a2b2"))
    
  ))
  
}

analyze_diversities <- function(phased_alleles_path, pipeline_path){
  
  type_seq <- sub(".*PHASE_([AB]+).*", "\\1", phased_alleles_path)

  AB_sequences <- suppressWarnings(toupper(read.fasta(gsub("phased", "unphased", phased_alleles_path), as.string = TRUE)))
  
  AB_alleles <- suppressWarnings(read.fasta(phased_alleles_path, as.string = TRUE))

  AB_alphas <- toupper(substr(AB_alleles, 1, 246))
  AB_betas <- toupper(substr(AB_alleles, 247, 519))
  
  AB_recombs <- get_recombs(AB_alphas, AB_betas)
  
  alpha_alleles_path <- sub("(.*PHASE_)AB", "\\1A", phased_alleles_path)
  beta_alleles_path <- sub("(.*PHASE_)AB", "\\1B", phased_alleles_path)
  
  alphas <- suppressWarnings(toupper(read.fasta(alpha_alleles_path, as.string = TRUE)))
  
  betas <- suppressWarnings(toupper(read.fasta(beta_alleles_path, as.string = TRUE)))
  
  alphabeta_recombs <- get_recombs(alphas, betas)
  
  individuals <- sub("[ab]$", "", names(AB_alleles))
  
  allele_comparisons_df <- data.frame(
    
    INDIVIDUAL = individuals,
    AB_ALPHAS = AB_alphas,
    ALPHAS = alphas,
    AB_BETAS = AB_betas,
    BETAS = betas
    
  )
  
  View(allele_comparisons_df)
  
  unique_individuals <- unique(individuals)
  
  cat("\nAlpha comparisons:\n")
  
  compare_alleles <- function(indv, ab_col, a_b_col, AB_sequences, A_B_sequences, init_index){
    
    indv_data <- allele_comparisons_df[allele_comparisons_df$INDIVIDUAL == indv, ]
    
    ab_data <- indv_data[[ab_col]]
    a_b_data <- indv_data[[a_b_col]]
    
    if(!(identical(sort(ab_data), sort(a_b_data)))){
      
      ab_chars_one <- strsplit(ab_data, "")[[1]]
      ab_chars_two <- strsplit(ab_data, "")[[2]]
      a_b_chars_one <- strsplit(a_b_data, "")[[1]]
      a_b_chars_two <- strsplit(a_b_data, "")[[2]]
      
      ab_seqs_chars <- strsplit(AB_sequences, "")[[1]]
      a_b_seqs_chars <- strsplit(A_B_sequences, "")[[1]]
      
      comparison_one <- which(ab_chars_one != a_b_chars_one)
      comparison_two <- which(ab_chars_one != a_b_chars_two)
      comparison_three <- which(ab_chars_two != a_b_chars_two)
      comparison_four <- which(ab_chars_two != a_b_chars_one)
      
      straight_mismatches <- unique(c(comparison_one, comparison_three))
      crossed_mismatches <- unique(c(comparison_two, comparison_four))
      
      #best_comparison <- list(straight_mismatches, crossed_mismatches)[[which.min(c(length(straight_mismatches), length(crossed_mismatches)))]]
      
      cat("\nINDV     BP    AB  A_B  AB_orig A_B_orig")
      cat("\n----------------------------------------\n")

      
      if(length(straight_mismatches) >= length(crossed_mismatches)){
        
        #comparison two is more parsimonious, use it here
        #hence, comp2 compares A1 with A4 and A2 with A3
        
        cat(paste(indv, comparison_two, ab_chars_one[comparison_two], a_b_chars_two[comparison_two], ab_seqs_chars[comparison_two + init_index], a_b_seqs_chars[comparison_two], sep="    "), sep="\n")
      
        
        #best_comparison_two <- comparison_four
        
        
        if(!identical(comparison_two, comparison_four) | 
           !identical(ab_chars_one[comparison_two], a_b_chars_one[comparison_four]) | 
           !identical(ab_chars_two[comparison_two], a_b_chars_two[comparison_four])){
          
          cat("/nSNPs likely in the other pair of alleles, here they are:/n")
          
          cat(paste(indv, comparison_four, ab_chars_two[comparison_four], a_b_chars_one[comparison_four], ab_seqs_chars[comparison_four + init_index], a_b_seqs_chars[comparison_four], sep="    "), sep="\n")
          
          
          
        }
        
        else{
          
          cat("The other alleles have no additional SNPs! Mismatches are the same and in the same locations!\n")
          
        }
        
        
      }
      
      else{ 
        
        #comparison one is more parsimonious, use it here
        #hence, comp1 compares A1 with A3 and A2 with A4
        
        cat(paste(indv, comparison_one, ab_chars_one[comparison_one], a_b_chars_one[comparison_one], ab_seqs_chars[comparison_one + init_index], a_b_seqs_chars[comparison_one], sep="    "), sep="\n")
        
        
        #best_comparison_two <- comparison_three
        
        if(!identical(comparison_one, comparison_three) | 
           !identical(ab_chars_one[comparison_one], a_b_chars_two[comparison_three]) | 
           !identical(ab_chars_two[comparison_one], a_b_chars_one[comparison_three])){
          
          cat("/nSNPs likely in the other pair of alleles, here they are:/n")
          
          cat(paste(indv, comparison_three, ab_chars_two[comparison_three], a_b_chars_two[comparison_three], ab_seqs_chars[comparison_three + init_index], a_b_seqs_chars[comparison_three], sep="    "), sep="\n")
          
          
        }
        
        else{
          
          cat("The other alleles have no additional SNPs! Mismatches are the same and in the same locations!\n")          
        }
        
        
      }

      
    }
    
  }
  
  sapply(unique_individuals, 
         compare_alleles, 
         ab_col = "AB_ALPHAS", 
         a_b_col = "ALPHAS", 
         AB_sequences = AB_sequences, 
         A_B_sequences = suppressWarnings(toupper(read.fasta(sub("phased", "unphased", alpha_alleles_path), as.string = TRUE))), 
         init_index=0)
  

  cat("\n\nBeta comparisons:\n\n")
  
  sapply(unique_individuals, 
         compare_alleles, 
         ab_col = "AB_BETAS", 
         a_b_col = "BETAS", 
         AB_sequences = AB_sequences, 
         A_B_sequences = suppressWarnings(toupper(read.fasta(sub("phased", "unphased", beta_alleles_path), as.string = TRUE))), 
         init_index=246)
  
  cat("\n")
  
  
}


main_genetic_checks <- function(phasepath, pipeline_path, p){
  
  if(!dir.exists(pipeline_path)) {
    dir.create(pipeline_path)
  }
  #creates folder for products of this pipeline if it doesn't already exist
  
  #summary stats of dataset --> dataframe/matrix --> save as CSV
  
  #problem is that some datasets were phased as all alphas or all betas, and some as AB together
  #best strategy is to PHASE A's and B's separately then concatenate them
  #what JLA did was PHASE AB's together, which loses some diversity due to PHASE pulling the alarm sooner
  
  comparative_out_data <- data.frame(rbind(
        #suppressWarnings(analyze_out(paste0(phasepath, "CP/PHASE_AB/seqphase.out"), pipeline_path, p)),
        suppressWarnings(analyze_out(paste0(phasepath, "p=0.5/NO_CLONES/PHASE_AB_NOCLONES/seqphase.out"), pipeline_path, p)),
        suppressWarnings(analyze_out(paste0(phasepath, "p=0.5/NO_CLONES/PHASE_A_NOCLONES/seqphase.out"), pipeline_path, p)),
        suppressWarnings(analyze_out(paste0(phasepath, "p=0.5/NO_CLONES/PHASE_B_NOCLONES/seqphase.out"), pipeline_path, p))
        #suppressWarnings(analyze_out(paste0(phasepath, "p=0.5/CLONES/PHASE_A_CLONES/seqphase.out"), pipeline_path, p)),
        #suppressWarnings(analyze_out(paste0(phasepath, "p=0.5/CLONES/PHASE_B_CLONES/seqphase.out"), pipeline_path, p)),
        #suppressWarnings(analyze_out(paste0(phasepath, "p=0.4/NO_CLONES/PHASE_AB_NOCLONES/seqphase.out"), pipeline_path, p)),
        #suppressWarnings(analyze_out(paste0(phasepath, "p=0.4/NO_CLONES/PHASE_A_NOCLONES/seqphase.out"), pipeline_path, p)),
        #suppressWarnings(analyze_out(paste0(phasepath, "p=0.4/NO_CLONES/PHASE_B_NOCLONES/seqphase.out"), pipeline_path, p)),
        #suppressWarnings(analyze_out(paste0(phasepath, "p=0.4/CLONES/PHASE_A_CLONES/seqphase.out"), pipeline_path, p)),
        #suppressWarnings(analyze_out(paste0(phasepath, "p=0.4/CLONES/PHASE_B_CLONES/seqphase.out"), pipeline_path, p))
        
        ))
  
  #View(comparative_out_data)
  
  #print(length(unique(toupper(read.fasta("/Users/mattracz/Projects/Wilson_Lab/Pipeline/Recombs.fasta", as.string=TRUE)))))
  #print(length(unique(toupper(read.fasta("/Users/mattracz/Projects/Wilson_Lab/Pipeline/recombs_postRECCO.fasta", as.string=TRUE)))))
  
  analyze_diversities(paste0(phasepath, "p=0.5/NO_CLONES/PHASE_AB_NOCLONES/phased.fasta"), pipeline_path)
  
  #after each out file processed, get number of unique recombinants from each AB outfile 
  #and from each Alpha + Beta outfiles in the same folder
  
  #get recombinants and their diversity (uniqueness, other diversity stat tests)
  
  #see if ABW's AB GS and my AB GS are the same 
  
  #check diversity of recombinants from AB phased versus AlphaBeta phased
  
  #get all recombinants from Alpha + Beta, then get number of uniques
  #also compare AB's alpha and beta to PHASEd alphas and betas
  #use from getUniqueRecombs function to make recombinants
  #compare diversity of recombinants 
  
  #if type_seq == AB --> make recombinants from alleles --> get number of unique recombinants 
  #if type_seq == A --> use for alphas and change A in pathname to B for betas --> make recombinants from alleles --> get number of unique recombinants 
  #compare alphas and betas to that PHASE's alpha and beta PHASEd alleles 
  #alphas and betas should be identical to AB's alphas and betas
  
}

main_genetic_checks("PHASED/", "PHASE_pipeline/", 0.7)
