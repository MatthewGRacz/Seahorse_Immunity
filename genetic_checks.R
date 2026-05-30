library("seqinr")

setwd("/Users/mattracz/Projects/Wilson_Lab")

analyze_out <- function(outpath, p){
  
  outfile <- readLines(outpath)
  
  type_seq <- sub(".*PHASE_([AB]+)_.*", "\\1", outpath)
  #grabs A, B, or AB in outpath name to classify allele as Alpha, Beta, or AlphaBeta
  
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
  
  phased_allele_freqs <- as.numeric(sub(".*\\s+([0-9]+\\.[0-9]+)\\s*$", "\\1", alleles_by_freq))
  #allele_numbers <- as.integer(sub("^\\s*([0-9]+)\\s+.*", "\\1", alleles_by_freq))
  
  #where says "positions of loci," get the loci positions
  
  allele_probs <- outfile[(grep("BEGIN PHASEPROBS", outfile, value=FALSE)+1):(grep("END PHASEPROBS", outfile, value=FALSE)-1)]
  
  allele_probs <- gsub("=", "1.00", allele_probs)
  
  allele_probs <- do.call(rbind, lapply(allele_probs, function(x){
    
    return(as.numeric(unlist(strsplit(x, "\\s+"))))
    
  }))
  
  loci_positions <- outfile[(grep("Positions of loci: ", outfile, value=FALSE))]
  loci_positions <- regmatches(loci_positions, gregexpr("[0-9]+", loci_positions))[[1]]
  
  colnames(allele_probs) <- loci_positions
  
  low_conf_loci <- I(apply(allele_probs, 1, function(x){
    
    loci <- names(x[x < p])
    
    if(length(loci)==0){
      
      return(0)
      
    }
      
    return(loci)
      
    
    
  }))
  
  input_alleles_data <- data.frame(
    
    INDIVIDUAL = sapply(alleles_by_indv, `[`, 1),
    NUM_LOW_CONFIDENCE = rowSums(allele_probs < p),
    LOW_CONFIDENCE_LOCI = low_conf_loci,
    NUM_LOCI = length(loci_positions),
    allele_probs,
    
    stringsAsFactors = FALSE,
    check.names = FALSE
    
  )
  
  View(input_alleles_data)
  
  phased_allele_seqs <- suppressWarnings(toupper(read.fasta(sub("seqphase.out", "phased.fasta", outpath), as.string=TRUE)))
  #get allelic sequences
  
  phased_allele_lengths <- nchar(phased_allele_seqs)
  
  phased_allele_numbers <- c(as.integer(sapply(alleles_by_indv, `[`, 2)), as.integer(sapply(alleles_by_indv, `[`, 3)))
  
  phased_allele_data <- data.frame(
    
    INDIVIDUAL = sapply(alleles_by_indv, `[`, 1),
    SEQUENCE = phased_allele_seqs,
    LENGTH = phased_allele_lengths,
    ALLELE = phased_allele_numbers,
    FREQUENCY = phased_allele_freqs
    
  )
  
  View(phased_allele_data)
  
  #get individual, allele1 sequence, its number, its frequency, then 
  
  #get number of flagged loci in outfile
  #get positions of loci in outfile
  #get number of individuals in outfile
  #get number of unique alleles in outfile
  #get p value used
  
  #"number of individuals"
  #"number of loci"
  #begin command line --> -p0.5
  #get number of probabilities below some threshold, name column "P>[threshold]" 
  #length(unique(allele_seqs))
  
  #return both as list, get [1] then [2] from it
  
  
}


main_genetic_checks <- function(phasepath, pipeline_path, p){
  
  #summary stats of dataset --> dataframe/matrix --> save as CSV
  
  #problem is that some datasets were phased as all alphas or all betas, and some as AB together
  #best strategy is to PHASE A's and B's separately then concatenate them
  #what JLA did was PHASE AB's together, which loses some diversity due to PHASE pulling the alarm sooner
  
  p <- 0.7
  
  suppressWarnings(analyze_out(paste0(phasepath, "p=0.5/NO_CLONES/PHASE_AB_NOCLONES/seqphase.out"), p))
  
  #after each out file processed, get number of unique recombinants from each AB outfile 
  #and from each Alpha + Beta outfiles in the same folder
  
  #get recombinants and their diversity (uniqueness, other diversity stat tests)
  
  #see if ABW's AB GS and my AB GS are the same 
  
  
}

main_genetic_checks("PHASED/", "PHASE_pipeline/", p)
