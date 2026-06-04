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

setwd("/Users/mattracz/Projects/Wilson_Lab")

get_outfile_analysis <- function(phasepath, p){
  
  outfile <- readLines(paste0(phasepath, "seqphase.out")) 
  
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
  
  return(list(input_alleles_data, phased_allele_freqs))
  
  
}

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

get_seq_analysis <- function(phasepath, alpha_phasepath, beta_phasepath, AB_out_freqs, A_out_freqs, B_out_freqs){
  
  phased_AB <- toupper(suppressWarnings(read.fasta(paste0(phasepath, "phased.fasta"), as.string=TRUE)))
  phased_A <- toupper(suppressWarnings(read.fasta(paste0(alpha_phasepath, "phased.fasta"), as.string=TRUE)))
  phased_B <- toupper(suppressWarnings(read.fasta(paste0(beta_phasepath, "phased.fasta"), as.string=TRUE)))
  
  unphased_AB <- suppressWarnings(toupper(unlist(read.fasta(paste0(phasepath, "unphased.fasta"), as.string=TRUE))))
  unphased_A <- suppressWarnings(unlist(read.fasta(paste0(alpha_phasepath, "unphased.fasta"), as.string=TRUE)))
  unphased_B <- suppressWarnings(unlist(read.fasta(paste0(beta_phasepath, "unphased.fasta"), as.string=TRUE)))
  
  individuals <- names(suppressWarnings(read.fasta(paste0(phasepath, "unphased.fasta"), as.string=TRUE)))
  
  return(data.frame(
    
    INDIVIDUAL = individuals,
    UNPHASED_AB = unphased_AB,
    AB_ALPHA = substr(phased_AB, 1, 246),
    AB_BETA = substr(phased_AB, 247, 519),
    AB_FREQUENCY = AB_out_freqs,
    UNPHASED_ALPHA = unphased_A,
    ALPHA = phased_A,
    ALPHA_FREQUENCY = A_out_freqs,
    UNPHASED_BETA = unphased_B,
    BETA = phased_B,
    BETA_FREQUENCY = B_out_freqs

    
  ))
  
  
}

main <- function(phasepath, p){
  
  alpha_phasepath <- sub("(.*PHASE_)AB", "\\1A", phasepath)
  beta_phasepath <- sub("(.*PHASE_)AB", "\\1B", phasepath)

  outfile <- readLines(paste0(phasepath, "seqphase.out")) 
  
  AB_out_analysis <- get_outfile_analysis(phasepath, p)
  AB_out_df <- AB_out_analysis[[1]]
  AB_out_freqs <- AB_out_analysis[[2]]
  #View(AB_out_df)

  A_out_analysis <- get_outfile_analysis(alpha_phasepath, p)
  A_out_df <- A_out_analysis[[1]]
  A_out_freqs <- A_out_analysis[[2]]
  #View(A_out_df)
  
  B_out_analysis <- get_outfile_analysis(beta_phasepath, p)
  B_out_df <- B_out_analysis[[1]]
  B_out_freqs <- B_out_analysis[[2]]
  #View(B_out_df)
  
  allele_seq_data_df <- get_seq_analysis(phasepath, alpha_phasepath, beta_phasepath, AB_out_freqs, A_out_freqs, B_out_freqs)
  View(allele_seq_data_df)
  
  mismatch_df <- get_allele_comparisons(allele_seq_data_df)
  View(mismatch_df)
  
}


main("PHASED/p=0.5/NO_CLONES/PHASE_AB_NOCLONES/", 0.7)








