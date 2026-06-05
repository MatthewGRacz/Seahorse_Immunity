# input: outfile
# get out loci etc for AB, A, and B outfiles --> ((A)(B))_out_analysis_df (input_alleles_data)
# send outpath to allele_assesser
# get frequency, allele_num for each allele
# changes outpath to phased --> gets A, B, AB_A, AB_B, AB alleles (make these all local vars)
# get original sequences for A, B, AB (change outpath to unphased) --> allele_seq_data_df (phased_allele_data + allele_comparison_df)
# send allele_seq_data_df to comparison function
# run all of the same comparison analyses using the data in the allele_seq_data_df dataset

library("seqinr")
library("waldo")
library("dplyr")
library("purrr")

setwd("/Users/mattracz/Projects/Wilson_Lab")

get_A_B_recombs <- function(alpha_seqs, beta_seqs){
  
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

get_AB_recombs <- function(AB_phasepath){
  
  
  indv_names <- unique(gsub("[ab]$", "", names(AB_phasepath)))
  #get names of individuals whose sequences you're looking at
  
  a1b1 <- toupper(paste0(substr(AB_phasepath[paste0(indv_names, "a")], 1, 246), substr(AB_phasepath[paste0(indv_names, "a")], 247, 519)))
  #a = 1st allele, b = 2nd allele, per individual
  #the snipped betas have 2 extra bps ahead and 1 extra bp behind, so cut them out
  #trims to functional reading frames
  a1b2 <- toupper(paste0(substr(AB_phasepath[paste0(indv_names, "a")], 1, 246), substr(AB_phasepath[paste0(indv_names, "b")], 247, 519)))
  a2b1 <- toupper(paste0(substr(AB_phasepath[paste0(indv_names, "b")], 1, 246), substr(AB_phasepath[paste0(indv_names, "a")], 247, 519)))
  a2b2 <- toupper(paste0(substr(AB_phasepath[paste0(indv_names, "b")], 1, 246), substr(AB_phasepath[paste0(indv_names, "b")], 247, 519)))
  
  return(c(
    
    setNames(a1b1, paste0(indv_names, "_a1b1")),
    setNames(a1b2, paste0(indv_names, "_a1b2")),
    setNames(a2b1, paste0(indv_names, "_a2b1")),
    setNames(a2b2, paste0(indv_names, "_a2b2"))
    
  ))
  
  
}

main <- function(phasepath, p){
  
  alpha_phasepath <- sub("(.*PHASE_)AB", "\\1A", phasepath)
  beta_phasepath <- sub("(.*PHASE_)AB", "\\1B", phasepath)
  
  AB_outpath <- paste0(phasepath, "seqphase.out")
  alpha_outpath <- paste0(alpha_phasepath, "seqphase.out")
  beta_outpath <- paste0(beta_phasepath, "seqphase.out")
  
  get_out_analysis <- function(outpath){
  
    outfile <- readLines(outpath) 
    #reads AB outfile
    
    allele_probs <- outfile[(grep("BEGIN PHASEPROBS", outfile, value=FALSE)+1):(grep("END PHASEPROBS", outfile, value=FALSE)-1)]
    #gets probabilities of the loci in the AB outfile
    
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
    
    alleles_by_indv <- outfile[(grep("BEGIN BESTPAIRS_SUMMARY", outfile, value=FALSE)+1):(grep("END BESTPAIRS_SUMMARY", outfile, value=FALSE)-1)]
    
    alleles_by_indv <- strsplit(gsub(":", ",", gsub("\\(|\\)", "", gsub(" ", "", alleles_by_indv))), ",")
    #[1] is name, [2] is allele1 number, [3] is allele2 number
    
    phased_allele_numbers <- c(rbind(as.integer(sapply(alleles_by_indv, `[`, 2)), as.integer(sapply(alleles_by_indv, `[`, 3))))
    
    phased_allele_freqs <- as.integer(sub(".*\\s+([0-9]+).*$", "\\1", outfile[(grep("BEGIN LIST_SUMMARY", outfile, value=FALSE)+1):(grep("END LIST_SUMMARY", outfile, value=FALSE)-1)]))[phased_allele_numbers]
    
    input_alleles_data <- data.frame(
      
      INDIVIDUAL = sapply(alleles_by_indv, `[`, 1),
      NUM_LOW_CONFIDENCE = rowSums(allele_probs < p),
      LOW_CONFIDENCE_LOCI = low_conf_loci,
      NUM_LOCI = length(loci_positions),
      allele_probs,
      
      stringsAsFactors = FALSE,
      check.names = FALSE
      
    )
    
    return(list(input_alleles_data, alleles_by_indv, phased_allele_freqs))
    
  }
  
  AB_out_data <- get_out_analysis(AB_outpath)
  AB_out_df <- AB_out_data[[1]]
  AB_allele_nums <- c(rbind(as.integer(sapply(AB_out_data[[2]], `[`, 2)), as.integer(sapply(AB_out_data[[2]], `[`, 3))))
  AB_freqs <- AB_out_data[[3]]
  
  alpha_out_data <- get_out_analysis(alpha_outpath)
  alpha_out_df <- alpha_out_data[[1]]
  alpha_allele_nums <- c(rbind(as.integer(sapply(alpha_out_data[[2]], `[`, 2)), as.integer(sapply(alpha_out_data[[2]], `[`, 3))))
  alpha_freqs <- alpha_out_data[[3]]
  
  beta_out_data <- get_out_analysis(beta_outpath)
  beta_out_df <- beta_out_data[[1]]
  beta_allele_nums <- c(rbind(as.integer(sapply(beta_out_data[[2]], `[`, 2)), as.integer(sapply(beta_out_data[[2]], `[`, 3))))
  beta_freqs <- beta_out_data[[3]]
  
  
  phased_AB <- suppressWarnings(toupper(read.fasta(paste0(phasepath, "phased.fasta"), as.string=TRUE)))
  phased_alpha <- suppressWarnings(toupper(read.fasta(paste0(alpha_phasepath, "phased.fasta"), as.string=TRUE)))
  phased_beta <- suppressWarnings(toupper(read.fasta(paste0(beta_phasepath, "phased.fasta"), as.string=TRUE)))
  
  unphased_AB <- suppressWarnings(toupper(read.fasta(paste0(phasepath, "unphased.fasta"), as.string=TRUE)))
  unphased_alpha <- suppressWarnings(toupper(read.fasta(paste0(alpha_phasepath, "unphased.fasta"), as.string=TRUE)))
  unphased_beta <- suppressWarnings(toupper(read.fasta(paste0(beta_phasepath, "unphased.fasta"), as.string=TRUE)))

  
  phased_df <- data.frame(
    
    INDIVIDUAL = rep(sapply(AB_out_data[[2]], `[`, 1), each = 2),
    
    AB_ALPHA = toupper(substr(phased_AB, 1, 246)),
    AB_BETA = toupper(substr(phased_AB, 247, 519)),
    AB_FREQUENCY = AB_freqs,
    
    ALPHA = phased_alpha,
    ALPHA_FREQUENCY = alpha_freqs,
    
    BETA = phased_beta,
    BETA_FREQUENCY = beta_freqs,
    
    stringsAsFactors = FALSE
    
    
  )
  
  unphased_df <- data.frame(
    
    INDIVIDUAL = sapply(AB_out_data[[2]], `[`, 1),
    UNPHASED_AB = unphased_AB,
    UNPHASED_ALPHA = unphased_alpha,
    UNPHASED_BETA = unphased_beta,
    
    stringsAsFactors = FALSE
    
  )
  

  
  allele_seq_df <<- merge(phased_df, unphased_df, by = "INDIVIDUAL", all.x = TRUE, sort = FALSE)
  
  allele_seq_df <- allele_seq_df[ , c(
    "INDIVIDUAL", "UNPHASED_AB", "AB_ALPHA", "AB_BETA", "AB_FREQUENCY", 
    "UNPHASED_ALPHA", "ALPHA", "ALPHA_FREQUENCY", 
    "UNPHASED_BETA", "BETA", "BETA_FREQUENCY" )]
  
  View(allele_seq_df)
  
  
  
  get_mismatches <- function(indv, allele_seq_df, ab_col, a_b_col, A_B_sequences, init_index, freq){
    
    indv_mismatches <- data.frame()
    
    indv_data <- allele_seq_df[allele_seq_df$INDIVIDUAL == indv, ]
    
    ab_data <- indv_data[[ab_col]]
    a_b_data <- indv_data[[a_b_col]]
    AB_sequences <- indv_data[["UNPHASED_AB"]]
    A_B_sequences <- indv_data[[A_B_sequences]]
    
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

      
      
       
      if(length(straight_mismatches) >= length(crossed_mismatches)){
        
        #comparison two is more parsimonious, use it here
        #hence, comp2 compares A1 with A4 and A2 with A3
        
        
        has_other_snps <- ((!identical(comparison_two, comparison_four)) ||
                         (!identical(ab_chars_one[comparison_two], a_b_chars_one[comparison_four])) || 
                         (!identical(ab_chars_two[comparison_two], a_b_chars_two[comparison_four])))
        
          
        
          indv_mismatches <- data.frame(
            
            INDIVIDUAL = indv,
            BP = comparison_two,
            AB = ab_chars_one[comparison_two],
            A_B = a_b_chars_two[comparison_two],
            AB_UNPHASED = ab_seqs_chars[comparison_two + init_index], 
            A_B_UNPHASED = a_b_seqs_chars[comparison_two],
            AB_FREQ = indv_data$AB_FREQUENCY[1],
            A_B_FREQ = indv_data[[freq]][2],
            OTHER_SNPS = has_other_snps,
            SECTION = a_b_col,
            
            stringsAsFactors = FALSE
            
          )
          
        #best_comparison_two <- comparison_four
        
        if(has_other_snps && !(is_empty(comparison_four)) && !(is_empty(comparison_two))){
          
            
            indv_mismatches <- data.frame(
              
              INDIVIDUAL = indv,
              BP = comparison_four,
              AB = ab_chars_two[comparison_four], 
              A_B = a_b_chars_one[comparison_four], 
              AB_UNPHASED = ab_seqs_chars[comparison_four + init_index], 
              A_B_UNPHASED = a_b_seqs_chars[comparison_four],
              AB_FREQ = indv_data$AB_FREQUENCY[2],
              A_B_FREQ = indv_data[[freq]][1],
              OTHER_SNPS = has_other_snps,
              SECTION = a_b_col,
              
              stringsAsFactors = FALSE
              
            )
          
        }
        

      }
      
      else{ 
        
        #comparison one is more parsimonious, use it here
        #hence, comp1 compares A1 with A3 and A2 with A4
        
        has_other_snps <- (!identical(comparison_one, comparison_three)) ||
                         (!identical(ab_chars_one[comparison_one], a_b_chars_two[comparison_three])) ||
                         (!identical(ab_chars_two[comparison_one], a_b_chars_one[comparison_three]))
        
          if(!(is_empty(comparison_one))){
          
          indv_mismatches <- data.frame(
            
            INDIVIDUAL = indv,
            BP = comparison_one,
            AB = ab_chars_one[comparison_one],
            A_B = a_b_chars_one[comparison_one],
            AB_UNPHASED = ab_seqs_chars[comparison_one + init_index],
            A_B_UNPHASED = a_b_seqs_chars[comparison_one],
            AB_FREQ = indv_data$AB_FREQUENCY[1],
            A_B_FREQ = indv_data[[freq]][1],
            OTHER_SNPS = has_other_snps,
            SECTION = a_b_col,
            
            stringsAsFactors = FALSE
            
            
          )
          }
        #best_comparison_two <- comparison_three
        
        if(has_other_snps && !(is_empty(comparison_three)) && !(is_empty(comparison_one))){
          
          
          
            indv_mismatches <- data.frame(
              
              INDIVIDUAL = indv,
              BP = comparison_three,
              AB = ab_chars_two[comparison_three], 
              A_B = a_b_chars_two[comparison_three], 
              AB_UNPHASED = ab_seqs_chars[comparison_three + init_index], 
              A_B_UNPHASED = a_b_seqs_chars[comparison_three], 
              AB_FREQ = indv_data$AB_FREQUENCY[2],
              A_B_FREQ = indv_data[[freq]][2],
              OTHER_SNPS = has_other_snps,
              SECTION = a_b_col,
              
              stringsAsFactors = FALSE
              
            )
          
        }
        
      }
      
      return(indv_mismatches)
      
    }
    
  }
  
  compare_alleles <- function(allele_seq_df, is_alpha){
    
    unique_individuals <- unique(allele_seq_df$INDIVIDUAL)
    
    alpha_mismatches <- data.frame()
    beta_mismatches <- data.frame()
    
    if(is_alpha){
      
      alpha_mismatches <- map_dfr(unique_individuals, function(indv) {
        
        get_mismatches(indv, allele_seq_df, "AB_ALPHA", "ALPHA", "UNPHASED_ALPHA", 0, "ALPHA_FREQUENCY")
        
      })
      
      return(alpha_mismatches)
      
      
    } else {
      

      
      beta_mismatches <- map_dfr(unique_individuals, function(indv) {
        
        get_mismatches(indv, allele_seq_df, "AB_BETA", "BETA", "UNPHASED_BETA", 246, "BETA_FREQUENCY")
      
        })
      
      return(beta_mismatches)
      
    }
    
    
  }
  
  mismatch_df <- data.frame(bind_rows(compare_alleles(allele_seq_df, is_alpha=TRUE), 
                           compare_alleles(allele_seq_df, is_alpha=FALSE)))
  
  View(mismatch_df)

  
}


main("PHASED/p=0.5/NO_CLONES/PHASE_AB_NOCLONES/", 0.7)








