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


setwd("/Users/mattracz/Projects/Wilson_Lab")

ABphasepath <- "PHASE_AB_2_17_2026/"
Aphasepath <- gsub("AB", "A", ABphasepath)
Bphasepath <- gsub("AB", "B", ABphasepath)

#----Clips AB sequence into Alpha/Beta:----

get_alpha <- function(seqs){
  
  cut_point <- 246
  
  num_sequences <- length(seqs)
  
  all_seqs <- list()
  all_names <- c()
  
  for(i in 1:num_sequences){ #for each unique allele
    
    #readseq <- read.fasta(seqs[i], as.string=TRUE)
    alpha <- toupper(substr(seqs[i], start=1, stop=cut_point))
    
    names(alpha) <- rep(names(seqs[i]), length(alpha))
    all_seqs <- c(all_seqs, alpha)           # Add sequences to the list
    all_names <- c(all_names, names(alpha))
    
    
  }
  
  write.fasta(sequences = all_seqs, names = all_names, file.out = "PHASE_A_2_16_2026/GS_167_A")
  
}
#gets alpha sections from GS 519bp sequences, puts them into FASTA file

get_alpha(read.fasta(paste0(ABphasepath, "GS_167_AB.fasta"), as.string=TRUE))

get_beta <- function(seqs){
  
  cut_point <- 247
  
  num_sequences <- length(seqs)
  
  all_seqs <- list()
  all_names <- c()
  
  for(i in 1:num_sequences){ #for each unique allele
    
    #readseq <- read.fasta(seqs[i], as.string=TRUE)
    beta <- toupper(substr(seqs[i], start=cut_point, stop=nchar(seqs[i])))
    
    names(beta) <- rep(names(seqs[i]), length(beta))
    all_seqs <- c(all_seqs, beta)           # Add sequences to the list
    all_names <- c(all_names, names(beta))
    
    
  }
  
  write.fasta(sequences = all_seqs, names = all_names, file.out = "PHASE_B_2_16_2026/GS_167_B")
  
}
#gets beta sections from GS 519bp sequences, puts them into FASTA file

get_beta(read.fasta(paste0(ABphasepath, "GS_167_AB.fasta"), as.string=TRUE))

#PHASE ALPHA WITH CLONES
#PHASE BETA WITH CLONES


#----Post-PHASE analaysis:----

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
#gets names of sequences in out file

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
    
    beginning = "BEGIN PHASEPROBS"
    
    probs_one <- c(as.numeric(unlist(strsplit(gsub("= ", "1.00 ", readout_one[sum(as.integer(grep(beginning, readout_one, value=FALSE)), as.integer(grep(indv_names[i], indv_names, value=FALSE)))]), " "))))
    #nth line down from "BEGIN PHASEPROBS" in the first out file, where n is individual n's position in the first fasta file
    
    
    for(k in 1:length(varpos_a_or_b)){ #1:1 correlation of positions of bps in sequence and probabilities in sequences
      
      if(probs_one[k] < prob){ 
        #is one of the alpha/beta out file loci mentioned the same position as a conflicting bp between two fasta sequences, or is its allelic position inferred with under 90% confidence?
        
        cat(paste0(indv_names[i], " at bp ", varpos_a_or_b[k], ": ", probs_one[k]))
        cat("\n")
        
      }
      
    }
    
    
  }
  
  
} 
#gets confidences of base calls in out file, flags any with less confidence than prob

getProbsOne(paste0(ABphasepath, "GS_167_AB.fasta"), 
            paste0(ABphasepath, "seqphase.out"),
            0.99)

#clones are phased due to being included in the total number of sequences in the input file, 
#which leads to N's and is redundant, so phased clones are excluded from AB phased alleles

getNoClonePath <- function(seqs){
  
  return (paste0(sub(".fasta", "", seqs), "_ClonesRemoved"))
  
} 
#returns the path name for the phased_NoClones file

removePhasedClones <- function(seqs){
  
  readseqs <- read.fasta(paste0(Aphasepath, "phased.fasta"), as.string=TRUE)
  
  uniqueseqnames <- names(readseqs)[!duplicated(names(readseqs))] #only names of the unique, non-duplicated seqs
  
  uniqueseqs <- lapply(readseqs[!duplicated(names(readseqs))], toupper) #only the seqs with the unique, non-duplicated names
  
  outpath <- getNoClonePath(seqs)
  
  write.fasta(sequences = uniqueseqs, names = uniqueseqnames, file.out = outpath)
  
  
} 
#removes clones from phased file, as they usually are redundant, not useful, and/or have N's 

suppressWarnings(removePhasedClones(Aphasepath))
suppressWarnings(removePhasedClones(Bphasepath))
suppressWarnings(removePhasedClones(ABphasepath))


hasNs <- function(seqs){
  
  readseqs <- read.fasta(getNoClonePath(paste0(seqs, "phased.fasta")), as.string=TRUE)
  
  totalmasks <- c()
  
  for(i in names(readseqs)){
      
    currentseq <- toupper(strsplit(readseqs[[i]], "")[[1]]) #all bases per sequence
    
    Nmask <- which(currentseq == "N") #mask where Ns appear in sequence
    
    totalmasks <- c(totalmasks, Nmask)
    
    if(!length(Nmask)==0){
      
      cat(paste0(i, " at bp(s): ", paste(Nmask, collapse = ", "))) #print bps where Ns appear
      cat("\n")
      
      
    }
    
  }
  
  if(length(totalmasks)==0){
    
    cat(paste0("No N's in ", paste0(Aphasepath, "phased.fasta"), "!")) #if no Ns in all masks
    cat("\n")
    return (TRUE)
    
  }
  
  return (FALSE) #if function gets here, didn't return TRUE, hence Nmask is not empty
  
  
} 
#checks if any non-clone alleles have N's in them

hasRetainedAlleles <- function(Aseqs, Bseqs, ABseqs){
  
  hasSameNumAlleles <- (((length(union((names(read.fasta(getNoClonePath(paste0(Aseqs, "phased.fasta")), as.string=TRUE))), (names(read.fasta(getNoClonePath(paste0(Bseqs, "phased.fasta")), as.string=TRUE)))))/2) - (length((names(read.fasta(getNoClonePath(paste0(ABseqs, "phased.fasta")), as.string=TRUE))))/2)) == 0)
  
  #checks if the same allele names are in A and B, and compares it to the number of allele names in AB; it subtracts the numbers of A_and_B alleles from AB alleles, and if it equals zero, they must have the same allele names and quantities, hence no alleles lost
  
  cat(paste0("Do AB and A and B have same number of alleles? ", hasSameNumAlleles))
  cat("\n")
  
  return(hasSameNumAlleles)

} 
#checks that no alleles were lost when clones removed


if( hasNs(Aphasepath) && hasNs(Aphasepath) && hasNs(Aphasepath) && hasRetainedAlleles(Aphasepath, Bphasepath, ABphasepath) ){
  
  cat("\n")
  cat("Files have no N's and no alleles lost, you may proceed")
  cat("\n")
  
} else {
  
  cat("\n")
  cat("Do NOT proceed")
  cat("Check previous messages for source of error")
  cat("\n")
  
}
#if no N's in alleles and no alleles lost, analysis may continue

getUniqueRecombs <- function(Aseqs, Bseqs, ABseqs){ #creates recombinants of AB sequences, gets number of unique AB recombs, gets number of unique combinations of Alpha and Beta alleles, compares them
  
  readAseqs <- read.fasta(getNoClonePath(paste0(Aseqs, "phased.fasta")), as.string=TRUE)
  readBseqs <- read.fasta(getNoClonePath(paste0(Bseqs, "phased.fasta")), as.string=TRUE)
  readABseqs <- read.fasta(getNoClonePath(paste0(ABseqs, "phased.fasta")), as.string=TRUE)
  
  cut_point = 246 #length of alpha sequence, bps after are beta
  
  A_B_recombs <- list()
  A_B_names <- c()
  
  AB_recombs <- list()
  AB_names <- c()
  
  for(i in 1:(length(readAseqs)/2)){ #lengths of A, B, and AB are same, so just using A for A and B, and AB for AB, for sake of clarity, but they're interchangeable
  #for recombs from Alpha and Beta sequences
    
    a1 <- toupper(readAseqs[[(2*i)-1]])[[1]]
    a2 <- toupper(readAseqs[[(2*i)]])[[1]]
    b1 <- toupper(readBseqs[[(2*i)-1]])[[1]]
    b2 <- toupper(readBseqs[[(2*i)]])[[1]]
    
    seqname <- gsub("[ab]$", "", names(readAseqs)[(2*i)-1]) #removes a/b from end of allele name
    
    a1b1 <- paste0(a1, b1)
    a2b1 <- paste0(a2, b1)
    a1b2 <- paste0(a1, b2)
    a2b2 <- paste0(a2, b2)
    
    A_B_recombs[[length(A_B_recombs) + 1]] <- a1b1
    A_B_recombs[[length(A_B_recombs) + 1]] <- a2b1
    A_B_recombs[[length(A_B_recombs) + 1]] <- a1b2
    A_B_recombs[[length(A_B_recombs) + 1]] <- a2b2
    
    A_B_names <- c(A_B_names, paste0(seqname, "_a1b1"))
    A_B_names <- c(A_B_names, paste0(seqname, "_a2b1"))
    A_B_names <- c(A_B_names, paste0(seqname, "_a1b2"))
    A_B_names <- c(A_B_names, paste0(seqname, "_a2b2"))
    
  } #Alpha alleles + Beta alleles recombinations
  
  A_B_unique <- unique(A_B_recombs)
  
  for(i in 1:(length(readABseqs)/2)){
    
    seq_one <- toupper(readABseqs[[(2*i)-1]])[[1]]
    seq_two <- toupper(readABseqs[[(2*i)]])[[1]]
    
    seqname <- gsub("[ab]$", "", names(readABseqs)[(2*i)-1])
    
    a1 <- substr(seq_one, 1, cut_point)
    a2 <- substr(seq_two, 1, cut_point)
    b1 <- substr(seq_one, cut_point+1, nchar(seq_one))
    b2 <- substr(seq_two, cut_point+1, nchar(seq_two))
    
    a1b1 <- paste0(a1, b1)
    a2b1 <- paste0(a2, b1)
    a1b2 <- paste0(a1, b2)
    a2b2 <- paste0(a2, b2)
    
    AB_recombs[[length(AB_recombs) + 1]] <- a1b1
    AB_recombs[[length(AB_recombs) + 1]] <- a2b1
    AB_recombs[[length(AB_recombs) + 1]] <- a1b2
    AB_recombs[[length(AB_recombs) + 1]] <- a2b2
    
    AB_names <- c(AB_names, paste0(seqname, "_a1b1"))
    AB_names <- c(AB_names, paste0(seqname, "_a2b1"))
    AB_names <- c(AB_names, paste0(seqname, "_a1b2"))
    AB_names <- c(AB_names, paste0(seqname, "_a2b2"))
    
    
  } #AB alleles recombinations
  
  AB_unique <- unique(AB_recombs)
  
  cat(paste0("Total number of unique recombinants that are in both AlphaBeta and AB: ", length(intersect(AB_unique, A_B_unique))))
  cat("\n")
  
  cat(paste0("Total number of unique recombinants exclusive to AlphaBeta: ", length(which(!(A_B_unique %in% AB_unique)))))
  cat("\n")
  
  cat(paste0("Total number of unique recombinants exclusive to AB: ", length(which(!(AB_unique %in% A_B_unique)))))
  cat("\n")
  
  cat(paste0("Total number of unique recombinants across both sets: ", length(union(AB_unique, A_B_unique)))) #union automatically removes duplicates, thus we get only the unique sequences that are in either list
  cat("\n")
  
  cat(paste0("Total number of unique recombinants in AlphaBeta: ", length(A_B_unique)))
  cat("\n")
  
  cat(paste0("Total number of unique recombinants in AB: ", length(AB_unique)))
  cat("\n")
  
  x=list("AB unique"=unlist(AB_unique), "AlphaBeta unique"=unlist(A_B_unique))
  
  maxval <- max(length(AB_unique), length(A_B_unique), length(union(AB_unique, A_B_unique)) )
  
  venndiag <- ggVennDiagram(x, 
                      label_alpha = 0, 
                      edge_size=0,
                      label = "both",
                      label_color = "white",
                      set_color = "white"
                      ) + 
    
    scale_fill_gradient2(low = "blue", mid = "hotpink", high = "red",
                         
                         limits= c(0, maxval),
                         midpoint= maxval/2,
                         name="Number of Alleles") +
    
    theme(plot.margin = margin(0, 20, 0, 50),
          plot.background = element_rect(fill = "navy", color = NA),
          legend.background = element_rect(fill = "skyblue"),
          legend.margin = margin(3, 3, 9, 3)) +
    
      coord_cartesian(clip = "off")
    
    print(venndiag)
  
  
} 
#gets number of unique AlphaBeta recombs, AB recombs, compares them, 
#and makes Venn Diagram of results

suppressMessages(suppressWarnings(getUniqueRecombs(Aphasepath, Bphasepath, ABphasepath)))



