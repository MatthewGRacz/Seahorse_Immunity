library("seqinr")

setwd("/Users/mattracz/Projects/Wilson_Lab")

analyze_out <- function(outpath){
  
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
  
  allele_freqs <- as.numeric(sub(".*\\s+([0-9]+\\.[0-9]+)\\s*$", "\\1", alleles_by_freq))
  allele_numbers <- as.integer(sub("^\\s*([0-9]+)\\s+.*", "\\1", alleles_by_freq))
  
  print(allele_numbers)
  
  phased_allele_seqs <- suppressWarnings(toupper(read.fasta(sub("seqphase.out", "phased.fasta", outpath), as.string=TRUE)))
  #get allelic sequences
  
  #print(phased_alleles)
  
  allele_lengths <- nchar(phased_allele_seqs)
  
  #where says "positions of loci," get the loci positions
  
  
  #View(allelic_info)
  
  
  
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


main_genetic_checks <- function(phasepath, pipeline_path){
  
  #summary stats of dataset --> dataframe/matrix --> save as CSV
  
  #problem is that some datasets were phased as all alphas or all betas, and some as AB together
  #best strategy is to PHASE A's and B's separately then concatenate them
  #what JLA did was PHASE AB's together, which loses some diversity due to PHASE pulling the alarm sooner
  
  analyze_out(paste0(phasepath, "p=0.5/NO_CLONES/PHASE_AB_NOCLONES/seqphase.out"))
  
  #after each out file processed, get number of unique recombinants from each AB outfile 
  #and from each Alpha + Beta outfiles in the same folder
  
  #get recombinants and their diversity (uniqueness, other diversity stat tests)
  
  #see if ABW's AB GS and my AB GS are the same 
  
  
}

main_genetic_checks("PHASED/", "PHASE_pipeline/")
