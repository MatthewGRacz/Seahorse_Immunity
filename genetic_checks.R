library("seqinr")

setwd("/Users/mattracz/Projects/Wilson_Lab")

analyze_out <- function(outpath, pipeline_path, p){
  
  outfile <- readLines(outpath)
  
  type_seq <- sub(".*PHASE_([AB]+)_.*", "\\1", outpath)
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
  
  #View(phased_allele_data)
  
  write.csv(phased_allele_data, paste0(pipeline_path, "phased_", type_seq, "_alleles_data.csv"))
  
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


main_genetic_checks <- function(phasepath, pipeline_path, p){
  
  #summary stats of dataset --> dataframe/matrix --> save as CSV
  
  #problem is that some datasets were phased as all alphas or all betas, and some as AB together
  #best strategy is to PHASE A's and B's separately then concatenate them
  #what JLA did was PHASE AB's together, which loses some diversity due to PHASE pulling the alarm sooner
  
  comparative_phase_data <- data.frame(rbind(suppressWarnings(analyze_out(paste0(phasepath, "p=0.5/NO_CLONES/PHASE_AB_NOCLONES/seqphase.out"), pipeline_path, p)),
        suppressWarnings(analyze_out(paste0(phasepath, "p=0.5/NO_CLONES/PHASE_A_NOCLONES/seqphase.out"), pipeline_path, p)),
        suppressWarnings(analyze_out(paste0(phasepath, "p=0.5/NO_CLONES/PHASE_B_NOCLONES/seqphase.out"), pipeline_path, p))))
  
  View(comparative_phase_data)
  
  print(length(unique(toupper(read.fasta("/Users/mattracz/Projects/Wilson_Lab/Pipeline/Recombs.fasta", as.string=TRUE)))))
  print(length(unique(toupper(read.fasta("/Users/mattracz/Projects/Wilson_Lab/Pipeline/recombs_postRECCO.fasta", as.string=TRUE)))))
  
  
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
