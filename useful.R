library("seqinr")
library("Biostrings")

setwd("/Users/mattracz/Projects/Wilson_Lab")

getOutNames <- function(readfile){
  
  beginning = "BEGIN BESTPAIRS_SUMMARY"
  ending = "END BESTPAIRS_SUMMARY"
  
  firstline <- as.integer(grep(beginning, readfile, value=FALSE)) +1
  lastline <- as.integer(grep(ending, readfile, value=FALSE)) -1
  
  outnames <- c()
  
  for(i in firstline:lastline){
    
    s <- readfile[i]
    outnames <- c(outnames, sub(":.*", "", s))
    
  }
  
  return (outnames)
  
}

getProbsOne <- function(fasta_one, out_one, prob){
  
  readfasta_one <- read.fasta(fasta_one)
  
  readout_one <- readLines(out_one)
  
  #find names in first out file, then second out file
  #then, get ith element (bc some names repeat) of indv_names (intersect of both, without CP) and see its position in first_names/second_names
  
  first_names <- unique(getOutNames(readout_one))#unique(gsub("[ab]$", "", names(readfasta_one)))
  #first_names <- first_names[!grepl("^CP", first_names)]
  #indv_names <- intersect(first_names, second_names)
  indv_names <- first_names[!grepl("^CP", first_names)]
  
  
  varpos_a_or_b <- sub("Positions of loci: ", "", readout_one[as.integer(grep("Positions of loci: ", readout_one, value=FALSE))])
  varpos_a_or_b <- as.integer(c(strsplit(varpos_a_or_b, " ")[[1]]))
  #being compared to other alpha/beta, so base positions remain same
  
  #variable positions, or loci where bps differ, among phased alleles in the first fasta file
  
  
  #<247 for alpha, >246 for beta; bps 1-246 of mhiiab sequences for alpha, bps 247-519 of mhiiab sequences for beta
  
  #variable positions, or loci where bps differ, among phased alleles in the second fasta file
  
  for(i in 1:length(indv_names)){#compare phased sequence of each shared individual to each other between fasta_one and fasta_two
    
    #bps <- c(which(as.character(readfasta_one[[paste0(indv_names[i], "a")]]) != as.character(readfasta_two[[paste0(indv_names[i], "a")]])))
    #bps <- bps[!is.na(bps)]
    #bp positions where sequences differ between both alpha/beta fasta files
    
    #if(length(bps) > 0){ #there are differences between the two fasta sequences for the same individual, individual n
    
    beginning = "BEGIN PHASEPROBS"
    
    probs_one <- c(as.numeric(unlist(strsplit(gsub("= ", "1.00 ", readout_one[sum(as.integer(grep(beginning, readout_one, value=FALSE)), as.integer(grep(indv_names[i], indv_names, value=FALSE)))]), " "))))
    #nth line down from "BEGIN PHASEPROBS" in the first out file, where n is individual n's position in the first fasta file
    
    #nth line down from "BEGIN PHASEPROBS" in the second out file, where n is individual n's position in the second fasta file
    
    #print(as.integer(grep(beginning, readout_one, value=FALSE)))
    #print(indv_names[i])
    #print(first_names)
    #print(as.integer(grep(indv_names[i], first_names, value=FALSE)))
    #print(as.integer(grep(indv_names[i], indv_names, value=FALSE)))
    #print(second_names)
    #print(as.integer(grep(indv_names[i], second_names, value=FALSE)))
    #print(indv_names)
    #print(as.integer(grep(indv_names[i], second_names, value=FALSE)))
    #print(probs_one)
    #print(probs_two)
    
    
    for(k in 1:length(varpos_a_or_b)){ #1:1 correlation of positions of bps in sequence and probabilities in sequences
      
      if(probs_one[k] < prob){ 
        #is one of the alpha/beta out file loci mentioned the same position as a conflicting bp between two fasta sequences, or is its allelic position inferred with under 90% confidence?
        
        print(paste0(indv_names[i], " at bp ", varpos_a_or_b[k], ": ", probs_one[k]))
        
      }
      
    }
    
    
    #}
    
  }
  
  
} #analyzes out file, flags base calls with confidence less than prob

getProbsOne("PHASE_AB_2_15_2026/GS_167_AB.fasta", 
            "PHASE_AB_2_15_2026/seqphase.out",
            0.99)

print(names(readDNAStringSet("PHASE_AB_2_15_2026/12_clones.fasta")))

