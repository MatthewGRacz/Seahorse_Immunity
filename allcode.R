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

#----JLA Clone and Captive Sequences Comparison Checker (in her format), also analyzes PHASE out file:----

getNumLoci <- function(file, after){
  
  return (as.integer(gsub(after, "", grep(after, readLines(file), value=TRUE))))
  
}

getListSummaryData <- function(file, beginning, ending){
  
    distinct_alleles <- c()
    seqs <- c()
    freqs <- c()
    one_off <- c()
    
    readfile <- readLines(file)  
  
    firstline <- as.integer(grep(beginning, readfile, value=FALSE)) + 1
    lastline <- as.integer(grep(ending, readfile, value=FALSE)) - 1
    
    
    for (i in firstline:lastline){ #for each 
      #seq_along(readfile)
      
      s <- readfile[i]
      s <- sub("^\\s+", "", s)
      n_allele <- as.integer(str_extract(s, "\\d+"))
      distinct_alleles <- c(distinct_alleles, n_allele)
      s <- sub(toString(n_allele), "", s)
      s <- sub("^\\s+", "", s)
      
      varpos_seq <- ""
      
      while(nchar(varpos_seq) < getNumLoci(file, "Number of Loci: ")){
        
        varpos_seq <- paste0(varpos_seq, str_extract(s, "\\d+"))
        s <- sub(str_extract(s, "\\d+"), "", s)
        s <- sub("^\\s+", "", s)
        
      }
      
      seqs <- c(seqs, varpos_seq)
      
      #print(varpos_seq)
      
      allele_freq <- as.integer(str_extract(s, "\\d+"))
      
      freqs <- c(freqs, allele_freq)
      if(allele_freq == 1){
        one_off <- c(one_off, TRUE)
      }
      else{
        one_off <- c(one_off, FALSE)
      }
      #print(allele_freq)
      
      
    }
    
    allele_data <- data.frame(distinct_alleles, seqs, freqs, one_off)
    return(allele_data)
    
  
}

MHIIAB_CP_data <- getListSummaryData("/Users/mattracz/Projects/2:5:21 Retests/MHIIAB-NOCLONES/CP/MHIIAB-CP.out", "BEGIN LIST_SUMMARY", "END LIST_SUMMARY")
View(MHIIAB_CP_data)

MHIIAB_data <- getListSummaryData("/Users/mattracz/Projects/2:5:21 Retests/MHIIAB-NOCLONES/NO CP/MHIIAB.out", "BEGIN LIST_SUMMARY", "END LIST_SUMMARY")
View(MHIIAB_data)


#BEGIN BESTPAIRS_SUMMARY get name of individual and two alleles, then see number of uniques in CP vs all else
#BEGIN PHASEPROBS get nth sequence, to that line, and get lowest and highest
#Get one-offs from BEST_PAIRS and freqs


getNumUniqueAlleles <- function(file){
  
  beginning = "BEGIN BESTPAIRS_SUMMARY"
  ending = "END BESTPAIRS_SUMMARY"
  
  readfile = readLines(file)
  
  nonCP_alleles <- c()
  CP_alleles <- c()
  
  firstline <- as.integer(grep(beginning, readfile, value=FALSE)) + 1
  lastline <- as.integer(grep(ending, readfile, value=FALSE)) - 1
  
  
  for (i in firstline:lastline){
    
    s <- readfile[i]
    indv_name <- str_split(s, ": ")[[1]][[1]]
    s <- sub(".*:", "", s)
    s <- sub("^\\s+", "", s)
    s <- sub("\\(", "", s)
    s <- sub("\\)", "", s)
    first_allele <- as.integer(str_split(s, "\\,")[[1]][[1]])
    second_allele <- as.integer(str_split(s, "\\,")[[1]][[2]])
    
    #print(first_allele)
    #print(second_allele)
    
    if(str_detect(indv_name, "CP")[1]){ #dataset has CP
      
      CP_alleles <- c(CP_alleles, first_allele, second_allele)
      
    }
    else{
      
      nonCP_alleles <- c(nonCP_alleles, first_allele, second_allele)
      
    }
    
  }

  return (c(length(unique(nonCP_alleles)),(length(setdiff(CP_alleles, nonCP_alleles))))) #number of unique alleles in main pop

}

getNumUniqueAlleles("/Users/mattracz/Projects/GS Phasing/ForST/GS.out")[[1]]

getNumUniqueAlleles("/Users/mattracz/Projects/2:5:21 Retests/MHIIA-CLONES/CP/MHIIA-CP-CLONES.out")[[1]]
getNumUniqueAlleles("/Users/mattracz/Projects/2:5:21 Retests/MHIIA-CLONES/CP/MHIIA-CP-CLONES.out")[[2]]
getNumUniqueAlleles("/Users/mattracz/Projects/2:5:21 Retests/MHIIA-CLONES/NO CP/MHIIA-CLONES.out")[[1]]
getNumUniqueAlleles("/Users/mattracz/Projects/2:5:21 Retests/MHIIA-CLONES/NO CP/MHIIA-CLONES.out")[[2]]

getNumUniqueAlleles("/Users/mattracz/Projects/2:5:21 Retests/MHIIB-CLONES/CP/MHIIB-CP-CLONES.out")[[1]]
getNumUniqueAlleles("/Users/mattracz/Projects/2:5:21 Retests/MHIIB-CLONES/CP/MHIIB-CP-CLONES.out")[[2]]
getNumUniqueAlleles("/Users/mattracz/Projects/2:5:21 Retests/MHIIB-CLONES/NO CP/MHIIB-CLONES.out")[[1]]
getNumUniqueAlleles("/Users/mattracz/Projects/2:5:21 Retests/MHIIB-CLONES/NO CP/MHIIB-CLONES.out")[[2]]

getNumUniqueAlleles("/Users/mattracz/Projects/2:5:21 Retests/MHIIAB-NOCLONES/CP/MHIIAB-CP.out")[[1]]
getNumUniqueAlleles("/Users/mattracz/Projects/2:5:21 Retests/MHIIAB-NOCLONES/CP/MHIIAB-CP.out")[[2]]
getNumUniqueAlleles("/Users/mattracz/Projects/2:5:21 Retests/MHIIAB-NOCLONES/NO CP/MHIIAB.out")[[1]]
getNumUniqueAlleles("/Users/mattracz/Projects/2:5:21 Retests/MHIIAB-NOCLONES/NO CP/MHIIAB.out")[[2]]

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

getOutNames(readLines("/Users/mattracz/Projects/2:5:21 Retests/MHIIA-CLONES/NO CP/MHIIA-CLONES.out"))

getOneOffs <- function(file, allele_data){
  
  beginning = "BEGIN BESTPAIRS_SUMMARY"
  ending = "END BESTPAIRS_SUMMARY"
  
  readfile = readLines(file)
  
  one_offs <- c()
  shared_alleles <- c()
  names <- c()
  
  firstline <- as.integer(grep(beginning, readfile, value=FALSE)) + 1
  lastline <- as.integer(grep(ending, readfile, value=FALSE)) - 1
  
  
  for (i in firstline:lastline){
    
    s <- readfile[i]
    indv_name <- str_split(s, ": ")[[1]][[1]]
    names <- c(names, indv_name)
    s <- sub(".*:", "", s)
    s <- sub("^\\s+", "", s)
    s <- sub("\\(", "", s)
    s <- sub("\\)", "", s)
    first_allele <- as.integer(str_split(s, "\\,")[[1]][[1]])
    second_allele <- as.integer(str_split(s, "\\,")[[1]][[2]])
    
    #if both one-offs, both into one_off and NULL for shared
    #else if both shared, both into shared_alleles and NULL for one_off
    #else if allele_one is one-off, into one_off and allele_two into shared_alleles
    #else allele_two into one_off and allele_one into shared_alleles
    
    if(allele_data[first_allele, "one_off"] && allele_data[second_allele, "one_off"]){
      one_offs <- c(one_offs, paste(first_allele, second_allele, sep=", "))
      shared_alleles <- c(shared_alleles, NA)
    }
    else if(!(allele_data[first_allele, "one_off"] || allele_data[second_allele, "one_off"])){
      shared_alleles <- c(shared_alleles, paste(first_allele, second_allele, sep=", "))
      one_offs <- c(one_offs, NA)
    }
    else if(allele_data[first_allele, "one_off"]){
      one_offs <- c(one_offs, first_allele)
      shared_alleles <- c(shared_alleles, second_allele)
    }
    else{
      one_offs <- c(one_offs, second_allele)
      shared_alleles <- c(shared_alleles, first_allele)
    }
    
  
  }
  
  return (data.frame(names, shared_alleles, one_offs))
  
}


onedata <- getOneOffs("/Users/mattracz/Projects/2:5:21 Retests/MHIIAB-NOCLONES/CP/MHIIAB-CP.out", MHIIAB_CP_data)
View(onedata)




suppressWarnings(getProbs("/Users/mattracz/Projects/2:5:21 Retests/TEST_ALPHA_CP",
                          "/Users/mattracz/Projects/2:5:21 Retests/MHIIA-CLONES/CP/MHIIA-CP-CLONES-PHASED.fasta",
                          "/Users/mattracz/Projects/2:5:21 Retests/MHIIA-CLONES/CP/MHIIA-CP-CLONES.out",
                          "/Users/mattracz/Projects/2:5:21 Retests/MHIIAB-NOCLONES/CP/MHIIAB-CP.out"
                          ))
#MHIIA vs ALPHA, with CP
suppressWarnings(getProbs("/Users/mattracz/Projects/2:5:21 Retests/TEST_ALPHA_NOCP",
                          "/Users/mattracz/Projects/2:5:21 Retests/MHIIA-CLONES/NO CP/MHIIA-CLONES-PHASED.fasta",
                          "/Users/mattracz/Projects/2:5:21 Retests/MHIIA-CLONES/NO CP/MHIIA-CLONES.out",
                          "/Users/mattracz/Projects/2:5:21 Retests/MHIIAB-NOCLONES/NO CP/MHIIAB.out"
                          ))
#MHIIA vs ALPHA, without CP ^
suppressWarnings(getProbs("/Users/mattracz/Projects/2:5:21 Retests/MHIIA-CLONES/NO CP/MHIIA-CLONES-PHASED.fasta",
                          "/Users/mattracz/Projects/2:5:21 Retests/MHIIA-CLONES/CP/MHIIA-CP-CLONES-PHASED.fasta",
                          "/Users/mattracz/Projects/2:5:21 Retests/MHIIA-CLONES/NO CP/MHIIA-CLONES.out",
                          "/Users/mattracz/Projects/2:5:21 Retests/MHIIA-CLONES/CP/MHIIA-CP-CLONES.out"))
#MHIIA vs MHIIA, no CP vs CP
suppressWarnings(getProbs("/Users/mattracz/Projects/2:5:21 Retests/MHIIB-CLONES/NO CP/MHIIB-CLONES-PHASED.fasta","/Users/mattracz/Projects/2:5:21 Retests/MHIIB-CLONES/CP/MHIIB-CP-CLONES-PHASED.fasta","/Users/mattracz/Projects/2:5:21 Retests/MHIIB-CLONES/NO CP/MHIIB-CLONES.out","/Users/mattracz/Projects/2:5:21 Retests/MHIIB-CLONES/CP/MHIIB-CP-CLONES.out"))
#MHIIB vs MHIIB, no CP vs CP
suppressWarnings(getProbs("/Users/mattracz/Projects/2:5:21 Retests/MHIIAB-NOCLONES/NO CP/MHIIAB-PHASED.fasta","/Users/mattracz/Projects/2:5:21 Retests/MHIIAB-NOCLONES/CP/MHIIAB-CP-PHASED.fasta","/Users/mattracz/Projects/2:5:21 Retests/MHIIAB-NOCLONES/NO CP/MHIIAB.out","/Users/mattracz/Projects/2:5:21 Retests/MHIIAB-NOCLONES/CP/MHIIAB-CP.out"))
#MHIIAB vs MHIIAB, no CP vs CP
suppressWarnings(getProbs("/Users/mattracz/Projects/2:5:21 Retests/TEST_BETA_NOCP",
                          "/Users/mattracz/Projects/2:5:21 Retests/MHIIB-CLONES/NO CP/MHIIB-CLONES-PHASED.fasta",
                          "/Users/mattracz/Projects/2:5:21 Retests/MHIIB-CLONES/NO CP/MHIIB-CLONES.out",
                          "/Users/mattracz/Projects/2:5:21 Retests/MHIIAB-NOCLONES/NO CP/MHIIAB.out"))
#MHIIB vs BETA, without CP #set alpha_length to 246
suppressWarnings(getProbs("/Users/mattracz/Projects/2:5:21 Retests/TEST_BETA_CP",
                          "/Users/mattracz/Projects/2:5:21 Retests/MHIIB-CLONES/CP/MHIIB-CP-CLONES-PHASED.fasta",
                          "/Users/mattracz/Projects/2:5:21 Retests/MHIIB-CLONES/CP/MHIIB-CP-CLONES.out",
                          "/Users/mattracz/Projects/2:5:21 Retests/MHIIAB-NOCLONES/CP/MHIIAB-CP.out"))
#MHIIB vs BETA, with CP #set alpha_length to 246



getProbs <- function(fasta_one, fasta_two, out_one, out_two){
  
  alpha_length = 0
  #set to 0 for alpha vs alpha, alpha vs mhiiab, mhiiab vs mhiiab, beta vs beta, but 246 for beta vs mhiiab
  
  readfasta_one <- read.fasta(fasta_one)
  readfasta_two <- read.fasta(fasta_two)
  
  readout_one <- readLines(out_one)
  readout_two <- readLines(out_two)
  
  #find names in first out file, then second out file
  #then, get ith element (bc some names repeat) of indv_names (intersect of both, without CP) and see its position in first_names/second_names
  
  first_names <- unique(getOutNames(readout_one))#unique(gsub("[ab]$", "", names(readfasta_one)))
  #first_names <- first_names[!grepl("^CP", first_names)]
  second_names <- unique(getOutNames(readout_two))#unique(gsub("[ab]$", "", names(readfasta_two)))
  #second_names <- second_names[!grepl("^CP", second_names)]
  indv_names <- intersect(first_names, second_names)
  indv_names <- indv_names[!grepl("^CP", indv_names)]
  
  
  varpos_a_or_b <- sub("Positions of loci: ", "", readout_one[as.integer(grep("Positions of loci: ", readout_one, value=FALSE))])
  varpos_a_or_b <- as.integer(c(strsplit(varpos_a_or_b, " ")[[1]]))
  #being compared to other alpha/beta, so base positions remain same

  #variable positions, or loci where bps differ, among phased alleles in the first fasta file
  
  varpos_ab_concat <- sub("Positions of loci: ", "", readout_two[as.integer(grep("Positions of loci: ", readout_two, value=FALSE))])
  varpos_ab_concat <- as.integer(c(strsplit(varpos_ab_concat, " ")[[1]]))
  varpos_ab_concat_filtered <- varpos_ab_concat[varpos_ab_concat>alpha_length] #for bps 247-519
  
  #<247 for alpha, >246 for beta; bps 1-246 of mhiiab sequences for alpha, bps 247-519 of mhiiab sequences for beta
  
  #variable positions, or loci where bps differ, among phased alleles in the second fasta file

  for(i in 1:length(indv_names)){#compare phased sequence of each shared individual to each other between fasta_one and fasta_two
    
    bps <- c(which(as.character(readfasta_one[[paste0(indv_names[i], "a")]]) != as.character(readfasta_two[[paste0(indv_names[i], "a")]])))
    bps <- bps[!is.na(bps)]
    #bp positions where sequences differ between both alpha/beta fasta files
    
    if(length(bps) > 0){ #there are differences between the two fasta sequences for the same individual, individual n
      
      beginning = "BEGIN PHASEPROBS"
      
      probs_one <- c(as.numeric(unlist(strsplit(gsub("= ", "1.00 ", readout_one[sum(as.integer(grep(beginning, readout_one, value=FALSE)), as.integer(grep(indv_names[i], first_names, value=FALSE)))]), " "))))
      #nth line down from "BEGIN PHASEPROBS" in the first out file, where n is individual n's position in the first fasta file
      
      
      
      
      probs_two <- c(as.numeric(unlist(strsplit(gsub("= ", "1.00 ", readout_two[sum(as.integer(grep(beginning, readout_two, value=FALSE)), as.integer(grep(indv_names[i], second_names, value=FALSE)))]), " "))))
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
      
      print("For MHIIA/B or without CP:")
      cat("\n")
      
      
      for(k in 1:length(varpos_a_or_b)){ #1:1 correlation of positions of bps in sequence and probabilities in sequences
        
        if(varpos_a_or_b[k] %in% bps || probs_one[k] < 0.50){ 
          #is one of the alpha/beta out file loci mentioned the same position as a conflicting bp between two fasta sequences, or is its allelic position inferred with under 90% confidence?
          
          print(paste0(indv_names[i], " at bp ", varpos_a_or_b[k], ": ", probs_one[k]))
          
        }
        
      }
      
      print(paste0("Maximum: ", max(probs_one)*100))
      print(paste0("Minimum: ", min(probs_one)*100))
      
      cat("\n")
      print("For MHIIAB or with CP:")
      cat("\n")
      
      
      
      beta_probs <- probs_two[as.integer(which(varpos_ab_concat == min(varpos_ab_concat_filtered))):length(probs_two)]
      #probabilities only for the beta half of the mhiiab sequences
      #if alpha/mhiiab, will just be a repeat of probs_two, as will j be a repeat of i
      
      
      for(k in varpos_ab_concat_filtered){#bps 247-519, correlate to positions in original varpos_ab_concat
        
        #print(i)
        
        j <- as.integer(which(varpos_ab_concat == k)) #get position of bp in probabilities sequence
        #j <- j[!is.na(j)]
        #print(j)
        
        
        if((k-alpha_length) %in% bps || probs_two[j] < 0.50){ 
          #is one of the mhiiab out file loci mentioned the same position as a conflicting bp between two fasta sequences, or is its corresponding allelic position inferred with under 90% confidence?
          #subtract 246 for beta because bp positions between two fasta sequences are 1-273, and beta bps in the mhiiab sequence are 247-519
          #so bp 1 for beta sequence is bp 247 for mhiiab sequence, etc
          
          print(paste0(indv_names[i], " at bp ", k-alpha_length, ": ", probs_two[j]))
          
        }
        
      }
      
      #print(paste0("Maximum: ", max(probs_two)*100))
      #print(paste0("Minimum: ", min(probs_two)*100))
      
      print(paste0("Maximum: ", max(beta_probs)*100))
      print(paste0("Minimum: ", min(beta_probs)*100))
      
      
      cat("\n")
      
      
    }
    
  }
  
  
} #USE THIS ONE



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
  
  
} #USE THIS ONE

getProbsOne("/Users/mattracz/Projects/Wilson_Lab/PHASE_AB_2_15_2026/GS_AB.fasta", 
            "/Users/mattracz/Projects/Wilson_Lab/PHASE_AB_2_15_2026/seqphase.out",
            0.99)

print(names(readDNAStringSet("/Users/mattracz/Projects/Wilson_Lab/PHASE_AB_2_15_2026/knowns.fasta")))


#----#RECCO, DataMonkey Analyses, Z-scores:----

#get functional opening reading frames for alleles (see frame where no stop codons in sequences, add N's/delete bases to accomodate)

recco <- function(alleles, reccodelim){
  
  read_alleles <- lapply(read.fasta(alleles, as.string = TRUE), toupper)
  
  recco_names <- read.delim(reccodelim)$Sequence
  
  kept_alleles <- read_alleles[!(names(read_alleles) %in% recco_names)]
  
  write.fasta(sequences = kept_alleles,names = names(kept_alleles),file.out = "/Users/mattracz/Projects/GS Phasing/Wilson_Lab/ForST/POST_RECCO_ALLELES")
  
  
}

recco("/Users/mattracz/Projects/Wilson_Lab/GS Phasing/ForST/RECOMB_STOPLESS_DNA_SEQS", "/Users/mattracz/Projects/GS Phasing/ForST/RECCO_results")

#DataMonkey -> identify which nucleotide sites are under positive selection
#get proteins to be supertyped

datamonkey <- function(meme, slac, fel){
  
  memeresults <- read_excel(meme)
  memesites <- na.omit(memeresults[memeresults$Class == "Diversifying", "Codon"]$Codon)
  
  felresults <- read_excel(fel)
  felsites <- na.omit(felresults[felresults$class == "Diversifying", "codon"]$codon)
  
  slacresults <- read.csv(slac)
  slacsites <- na.omit(slacresults[slacresults$`P..dN.dS...1.` < 0.05, ]$Site)
  
  allsites <- c(memesites, felsites, slacsites)
  codons <- as.numeric(names(table(allsites)[table(allsites) > 1]))
  
  usedalleles <- read.fasta("/Users/mattracz/Projects/GS Phasing/ForST/POST_RECCO_ALLELES")
  usedallele_names <- names(usedalleles)
  
  proteins <- lapply(usedalleles, translate)
  
  aminoacids <- list()
  
  for(n in usedallele_names){
    
    aas <- proteins[[n]][codons]
    aas <- paste(aas, collapse="")
    aminoacids[[length(aminoacids)+1]] <- aas
    
  }
  
  write.fasta(sequences=aminoacids, names=usedallele_names, file.out="/Users/mattracz/Projects/GS Phasing/ForST/AMINOACIDS")
  
}

datamonkey("/Users/mattracz/Projects/MEME_results.xlsx", "/Users/mattracz/Projects/SLAC_results.csv", "/Users/mattracz/Projects/FEL_results.xlsx")

#split each protein into individual amino acids, get the 5 Z values for each amino acid
zanalysis <- function(proteins){
  
  readproteins <- read.fasta("/Users/mattracz/Projects/GS Phasing/ForST/AMINOACIDS", as.string=TRUE)
  
  prots <- c()
  zlogs <- c()
  
  for(n in names(readproteins)){
    
    charprot <- strsplit(readproteins[[n]], "")[[1]]
    
    prots[[n]] <- (charprot)
    
    aazs <- lapply(charprot, zScales)
    
    zlogs[[n]] <- unlist(aazs)
    
  }
  
  zscdata <- data.frame(Allele=names(zlogs),
                        do.call(rbind, zlogs),
                        row.names=NULL
  )
  
  View(zscdata)
  
  return(zscdata)
  
  
  
  
  #JLA
  
  #Code2split <- data.frame(matrix(unlist(prots), nrow=length(prots), byrow=T))
  #output<-sapply(Code2split,zScales)
  #output2<-as.data.frame(do.call(rbind, output))
  #tmp <- data.frame(
  #  output2,
  #  ind=rep(1:nrow(output), nrow(output2)/nrow(output)),
  #  aa=rep(1:(nrow(output2)/nrow(output)), each=nrow(output))   
  #)
  #test<-reshape(tmp, idvar = "ind", timevar = "aa", direction = "wide")
  #row.names(test)<-paste(names(prots))
  #Testclusters <- find.clusters(test, max.n.clust=40)
  
  #View(test)
  
  #JLA
  
  
  
}

zscdata <- zanalysis("/Users/mattracz/Projects/GS Phasing/ForST/AMINOACIDS")

#DAPC supertyping (clustering) analysis for proteins

zscores_unlabeled <- zscdata[, grepl("^Z", names(zscdata))]

View(zscores_unlabeled)


#----Allele Subtraction:----

origDNA <- "ATGGCTACCTGTKKCAYRYGGACGKTGGTTGTGTGTTCAACTCGAGTGACCYGAATGACATCGAGTACTTCCAGATTTACAACTACAACAAACTGAAGCTTTTCCGCTTCAGCAGCACTTTGGATAAGTACGTCGGCTACACMGAGTTTGGCATCAAGCAGGCTACCGCCTTCAACAACRACAAAGACWWCATCGCCGACGYCAGAGCCAKGAAAGAAYACCTTTGTTTAAACAATATTAAGMWTGACTACGAAAGTGCGCTCACCAAGTCAG"

origDNA <- toupper(origDNA)
origDNA <- strsplit(origDNA, "")[[1]]

first_allele <- "ATGGCTACCTGTTTCATACGGACGGTGGTTGTGTGTTCAACTCGAGTGACCCGAATGACATCGAGTACTTCCAGATTTACAACTACAACAAACTGAAGCTTTTCCGCTTCAGCAGCACTTTGGATAAGTACGTCGGCTACACAGAGTTTGGCATCAAGCAGGCTACCGCCTTCAACARCGACAAAGACATCATCGCCGACGTCAGAGCCAGGAAAGAATACCTTTGTTTAAACAATATTAAGCTTGACTACGAAAGTGCGCTCACCAAGTCAG"

first_allele <- toupper(first_allele)
first_allele <- strsplit(first_allele, "")[[1]]

second_allele <- origDNA #if one-off error in allele, we know the answer in the original DNA sequence

varcodes <- list("M", "K", "Y", "R", "S", "W")

origbps <- c(origDNA[(which(origDNA != first_allele))])
print(origbps)
allelebps <- c(first_allele[(which(first_allele != origDNA))])
print(allelebps)
diff_indices <- (which(origDNA != first_allele))
print(diff_indices)


for(i in 1:length(diff_indices)){
  
  if((allelebps[i] %in% varcodes) && !(origbps[i] %in% varcodes)){ #heterozygotic label in allele, but not in original sequence
    
    print(paste0("FLAGGED: at bp ", diff_indices[i]))
    print(paste0("in the DNA sequence: ", origbps[i]))
    print(paste0("in the allele sequence: ", allelebps[i]))
    cat("\n")
  }
  
  else if((origbps[i] %in% varcodes) && !(allelebps[i] %in% varcodes)){ #heterozygotic label in original sequence, and not in allele
    
    if(origbps[i] == "K" && allelebps[i] == "T"){
      
      second_allele[diff_indices[i]] <- "G"
      
    }
    
    else if(origbps[i] == "K" && allelebps[i] == "G"){
      
      second_allele[diff_indices[i]] <- "T"
      
    }
    
    else if(origbps[i] == "Y" && allelebps[i] == "T"){
      
      second_allele[diff_indices[i]] <- "C"
      
    }
    
    else if(origbps[i] == "Y" && allelebps[i] == "C"){
      
      second_allele[diff_indices[i]] <- "T"
      
    }
    
    else if(origbps[i] == "S" && allelebps[i] == "G"){
      
      second_allele[diff_indices[i]] <- "C"
      
    }
    
    else if(origbps[i] == "S" && allelebps[i] == "C"){
      
      second_allele[diff_indices[i]] <- "G"
      
    }
    
    else if(origbps[i] == "W" && allelebps[i] == "T"){
      
      second_allele[diff_indices[i]] <- "A"
      
    }
    
    else if(origbps[i] == "W" && allelebps[i] == "A"){
      
      second_allele[diff_indices[i]] <- "T"
      
    }
    
    else if(origbps[i] == "M" && allelebps[i] == "A"){
      
      second_allele[diff_indices[i]] <- "C"
      
    }
    
    else if(origbps[i] == "M" && allelebps[i] == "C"){
      
      second_allele[diff_indices[i]] <- "A"
      
    }
    
    else if(origbps[i] == "R" && allelebps[i] == "A"){
      
      second_allele[diff_indices[i]] <- "G"
      
    }
    
    else if(origbps[i] == "R" && allelebps[i] == "G"){
      
      second_allele[diff_indices[i]] <- "A"
      
    }
    
    
  }
  
  
  
  
}

print(paste0("ALLELE 1: ", paste(first_allele, collapse="")))
print(paste0("ALLELE 2 (VIA SUBTRACTION): ", paste(second_allele, collapse="")))



#----Microbe Database and RNA Mapping (WIP):----

Sys.setenv(PATH = paste("/Users/mattracz/Projects/Microbial_real/ncbi-blast-2.16.0+/bin", Sys.getenv("PATH"), sep=":"))
#allows BLAST+ package to run
setwd("/Users/mattracz/Projects/Microbial_real")
#system("which blastn")
#system("which blastp")

blast_setup <- blast(db="/Users/mattracz/Projects/Microbial_real/BLAST+databases/16S_ribosomal_RNA/16S_ribosomal_RNA", remote=FALSE, type="blastn")
#sets up blast database for 16S RNA

get_order <- function(){
  
  #files in AA_seqs that are not directories (actual files)
  files <- list.files("/Users/mattracz/Projects/Microbial_real/RNA_seqs/", full.names=TRUE)[file.info(list.files("/Users/mattracz/Projects/Microbial_real/RNA_seqs/", full.names=TRUE))$isdir == FALSE]  # Filter out directories
  
  for(s in files){
    #for each AA sequence just created and stored in AA_seqs
    current_seq <- read.fasta(s) #AA sequence used
    current_seq <- paste(current_seq[[1]], collapse = "") #get AA sequence as an element
    
    RNAseq <- RNAStringSet(x=current_seq, start=1, end=nchar(current_seq), width=NA, use.names=TRUE)
    #make the AA sequence an AAStringSet so that BLAST can operate on it
    
    blast_data <- predict(blast_setup, RNAseq) #gets data of protein
    #print(blast_data) #print results of blast search on protein
    
    if(nrow(blast_data)==0){ #if no match
      print("nothing")
      individual_file_name <- strsplit(basename(s), "/Users/mattracz/Projects/Microbial_real/RNA_seqs/")[[1]]
      file.rename(from=s, to=paste("/Users/mattracz/Projects/Microbial_real/RNA_seqs/rejected/", individual_file_name,  sep=""))
      next
    }
    
    
    taxsearch <- entrez_search(db="nuccore", term=blast_data$sseqid[1])
    tax_id <- entrez_summary(db="nuccore", id=taxsearch$ids[1])$taxid
    #print(tax_id)#print that taxonomic ID
    
    taxonomic_info <- entrez_fetch(db="taxonomy", id=tax_id, rettype="xml", parsed=TRUE)
    #then fetch the full taxonomic lineage from that record
    
    #print(taxonomic_info) #View full taxonomic lineage
    
    order <- xpathSApply(taxonomic_info, "//Taxon[Rank='order']/ScientificName", xmlValue)
    #get the family of its taxonomic lineage and convert it to text content (not necessarily a string, might be a character code)
    
    #print(family) #print the family's name
    
    print(order)
    
    individual_file_name <- strsplit(basename(s), "/Users/mattracz/Projects/Microbial_real/RNA_seqs/")[[1]]
    file.rename(from=s, to=paste("/Users/mattracz/Projects/Microbial_real/RNA_seqs/kept/", individual_file_name,  sep=""))
    
    
  } 
  
}

get_order()





#----Microbe Database and Protein Mapping (WIP):----

Sys.setenv(PATH = paste("/Users/mattracz/Projects/Microbial_real/ncbi-blast-2.16.0+/bin", Sys.getenv("PATH"), sep=":"))
#allows BLAST+ package to run

#system("which blastn")
#system("which blastp")

blast_setup <- blast(db = "/Users/mattracz/Projects/Microbial_real/BLAST+databases/nr/nr.000", remote = FALSE, type = "blastp")
#sets up blast database for 16S RNA

get_order <- function(){
  
  #files in AA_seqs that are not directories (actual files)
  files <- list.files("/Users/mattracz/Projects/Microbial_real/AA_seqs/", full.names=TRUE)[file.info(list.files("/Users/mattracz/Projects/Microbial_real/AA_seqs/", full.names=TRUE))$isdir == FALSE]  # Filter out directories
  
  for(s in files){
    #for each AA sequence just created and stored in AA_seqs
    current_seq <- read.fasta(s) #AA sequence used
    current_seq <- paste(current_seq[[1]], collapse = "") #get AA sequence as an element
    
    AAseq <- AAStringSet(x=current_seq, start=1, end=nchar(current_seq), width=NA, use.names=TRUE)
    #make the AA sequence an AAStringSet so that BLAST can operate on it
    
    blast_data <- predict(blast_setup, AAseq) #gets data of protein
    #print(blast_data) #print results of blast search on protein
    
    if(nrow(blast_data)==0){ #if no match
      print("nothing")
      individual_file_name <- strsplit(basename(s), "/Users/mattracz/Projects/Microbial_real/AA_seqs/")[[1]]
      file.rename(from=s, to=paste("/Users/mattracz/Projects/Microbial_real/AA_seqs/rejected/", individual_file_name,  sep=""))
      next
    }
    
    taxsearch <- entrez_search(db="protein", term=blast_data$sseqid[1])
    
    if(length(taxsearch$ids)==0 || is.null(taxsearch$ids[1])){
      print("nothing")
      individual_file_name <- strsplit(basename(s), "/Users/mattracz/Projects/Microbial_real/AA_seqs/")[[1]]
      file.rename(from=s, to=paste("/Users/mattracz/Projects/Microbial_real/AA_seqs/rejected/", individual_file_name,  sep=""))
      next
    }
    
    tax_id <- entrez_summary(db="protein", id=taxsearch$ids[1])$taxid
    
    
    #print(tax_id)#print that taxonomic ID
    
    taxonomic_info <- entrez_fetch(db="taxonomy", id=tax_id, rettype="xml", parsed=TRUE)
    #then fetch the full taxonomic lineage from that record
    
    #print(taxonomic_info) #View full taxonomic lineage
    
    order <- xpathSApply(taxonomic_info, "//Taxon[Rank='order']/ScientificName", xmlValue)
    #get the family of its taxonomic lineage and convert it to text content (not necessarily a string, might be a character code)
    
    #print(family) #print the family's name
    
    print(order)
    
    individual_file_name <- strsplit(basename(s), "/Users/mattracz/Projects/Microbial_real/AA_seqs/")[[1]]
    file.rename(from=s, to=paste("/Users/mattracz/Projects/Microbial_real/AA_seqs/kept/", individual_file_name,  sep=""))
    
  }
  
}

get_order()



#----General Use Protein, adds FASTA sequences, and AB clipper:----

all_seqs <- list()
all_names <- c()

combine_DNA_allele_seqs <- function(seqs){
  for (f in list.files(seqs, full.names=TRUE)){
    #seqs <- read.fasta(f, as.string = TRUE)
    seqname <- basename(f)
    #print(seqname)
    aa <- toupper(read.fasta(f, as.string = TRUE))
    names(aa) <- rep(seqname, length(aa))
    all_seqs <- c(all_seqs, aa)           # Add sequences to the list
    all_names <- c(all_names, names(aa))
  }
  
  write.fasta(sequences = all_seqs, names = all_names, file.out = "/Users/mattracz/Projects/Microbial_real/SINGLES_DNA_SEQS")
}

combine_RNA_allele_seqs <- function(seqs){
  for (f in list.files(seqs, full.names=TRUE)){
    #seqs <- read.fasta(f, as.string = TRUE)
    seqname <- basename(f)
    #print(seqname)
    aa <- toupper(read.fasta(f, as.string = TRUE))
    names(aa) <- rep(seqname, length(aa))
    all_seqs <- c(all_seqs, aa)           # Add sequences to the list
    all_names <- c(all_names, names(aa))
  }
  
  write.fasta(sequences = all_seqs, names = all_names, file.out = "/Users/mattracz/Projects/Microbial_real/FINALSEQS_RNA")
}

combine_AA_allele_seqs <- function(seqs){
  for (f in list.files(seqs, full.names=TRUE)){
    #seqs <- read.fasta(f, as.string = TRUE)
    seqname <- basename(f)
    #print(seqname)
    aa <- toupper(read.fasta(f, as.string = TRUE))
    names(aa) <- rep(seqname, length(aa))
    all_seqs <- c(all_seqs, aa)           # Add sequences to the list
    all_names <- c(all_names, names(aa))
  }
  
  write.fasta(sequences = all_seqs, names = all_names, file.out = "/Users/mattracz/Projects/Microbial_real/FINALSEQS_AA")
}

combine_DNA_allele_seqs("/Users/mattracz/Projects/Microbial_real/DNA_seqs/singles")
#combine_RNA_allele_seqs("/Users/mattracz/Projects/Microbial_real/RNA_seqs/kept")
#combine_AA_allele_seqs("/Users/mattracz/Projects/Microbial_real/AA_seqs/kept")

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
  
  write.fasta(sequences = all_seqs, names = all_names, file.out = "/Users/mattracz/Projects/TEST_ALPHA_CP")
  
}

get_alpha(read.fasta("/Users/mattracz/Projects/MHIIAB-CP-PHASED.fasta", as.string=TRUE))


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
  
  write.fasta(sequences = all_seqs, names = all_names, file.out = "/Users/mattracz/Projects/TEST_BETA_CP")
  
}

get_beta(read.fasta("/Users/mattracz/Projects/MHIIAB-CP-PHASED.fasta", as.string=TRUE))


#----Makes AB recombinations, converts stop codon bases to Ns, converts to proteins, checks if Sygnathid, keeps only codons identified as under positive selection (via DataMonkey), gets thwir Z-scores and runs DAPC analyses (WIP):----

#----Generates Random RNA, DNA, and Protein FASTA Sequences:----

nfiles = 10

rna_generate <-function(rand_sequence, seqlength){
  for(i in 1:seqlength) {
    
    n=sample(1:4, 1, replace=TRUE)
    
    if(n ==1) {
      rand_sequence <- c(rand_sequence, "A")  # Add "A" to the sequence
    } else if(n ==2) {
      rand_sequence <- c(rand_sequence, "G")  # Add "G" to the sequence
    } else if(n ==3) {
      rand_sequence <- c(rand_sequence, "U")  # Add "U" to the sequence
    } else if(n ==4) {
      rand_sequence <- c(rand_sequence, "C")  # Add "C" to the sequence
    }
  }
  return(rand_sequence)
}

dna_generate <-function(rand_sequence, seqlength){
  for(i in 1:seqlength) {
    
    n=sample(1:4, 1, replace=TRUE)
    
    if(n ==1) {
      rand_sequence <- c(rand_sequence, "A")  # Add "A" to the sequence
    } else if(n ==2) {
      rand_sequence <- c(rand_sequence, "G")  # Add "G" to the sequence
    } else if(n ==3) {
      rand_sequence <- c(rand_sequence, "T")  # Add "T" to the sequence
    } else if(n ==4) {
      rand_sequence <- c(rand_sequence, "C")  # Add "C" to the sequence
    }
  }
  
  return(rand_sequence)
}

aa_generate <-function(rand_sequence, seqlength){
  for(i in 1:seqlength) {
    
    n=sample(1:20, 1, replace=TRUE)
    
    if(n ==1) {
      rand_sequence <- c(rand_sequence, "A")  # Add "A" to the sequence
    } else if(n ==2) {
      rand_sequence <- c(rand_sequence, "C")  # Add "C" to the sequence
    } else if(n ==3) {
      rand_sequence <- c(rand_sequence, "D")  # Add "D" to the sequence
    } else if(n ==4) {
      rand_sequence <- c(rand_sequence, "E")  # Add "E" to the sequence
    }else if(n ==5) {
      rand_sequence <- c(rand_sequence, "F")  # Add "F" to the sequence
    } else if(n ==6) {
      rand_sequence <- c(rand_sequence, "G")  # Add "G" to the sequence
    } else if(n ==7) {
      rand_sequence <- c(rand_sequence, "H")  # Add "H" to the sequence
    }else if(n ==8) {
      rand_sequence <- c(rand_sequence, "I")  # Add "I" to the sequence
    } else if(n ==9) {
      rand_sequence <- c(rand_sequence, "K")  # Add "K" to the sequence
    } else if(n ==10) {
      rand_sequence <- c(rand_sequence, "L")  # Add "L" to the sequence
    }else if(n ==11) {
      rand_sequence <- c(rand_sequence, "M")  # Add "M" to the sequence
    } else if(n ==12) {
      rand_sequence <- c(rand_sequence, "N")  # Add "N" to the sequence
    } else if(n ==13) {
      rand_sequence <- c(rand_sequence, "P")  # Add "P" to the sequence
    }else if(n ==14) {
      rand_sequence <- c(rand_sequence, "Q")  # Add "Q" to the sequence
    } else if(n ==15) {
      rand_sequence <- c(rand_sequence, "R")  # Add "R" to the sequence
    } else if(n ==16) {
      rand_sequence <- c(rand_sequence, "S")  # Add "S" to the sequence
    }else if(n ==17) {
      rand_sequence <- c(rand_sequence, "T")  # Add "T" to the sequence
    } else if(n ==18) {
      rand_sequence <- c(rand_sequence, "V")  # Add "V" to the sequence
    } else if(n ==19) {
      rand_sequence <- c(rand_sequence, "W")  # Add "W" to the sequence
    }else if(n ==20) {
      rand_sequence <- c(rand_sequence, "Y")  # Add "Y" to the sequence
    }
  }
  
  return(rand_sequence)
}

run <- function(seq_type, seqlength){
  
  for(filenum in 1:nfiles){
    rand_sequence <- character()  # Create an empty character vector for the RNA sequence
    seqlength = seqlength
    
    filepath <- ""
    
    # Convert the character vector to a single string
    rand_string <- paste(rand_sequence, collapse = "")
    
    if(seq_type == 1){
      filepath <- paste("/Users/mattracz/Projects/Microbial_real/RNA_seqs/", "random_gen_sequence", filenum, ".fasta", sep="")
      rand_sequence <- rna_generate(rand_sequence, seqlength)
      rand_string <- paste(rand_sequence, collapse = "")
      
    } else if(seq_type == 2){
      filepath <- paste("/Users/mattracz/Projects/Microbial_real/DNA_seqs/", "random_gen_sequence", filenum, ".fasta", sep="")
      rand_sequence <- dna_generate(rand_sequence, seqlength)
      rand_string <- paste(rand_sequence, collapse = "")
      
    } else if(seq_type == 0){
      filepath <- paste("/Users/mattracz/Projects/Microbial_real/AA_seqs/", "random_gen_sequence", filenum, ".fasta", sep="")
      rand_sequence <- aa_generate(rand_sequence, seqlength)
      rand_string <- paste(rand_sequence, collapse = "")
      
    } 
    #print(filepath)
    
    # Write the RNA sequence to a FASTA file
    write.fasta(sequence=rand_string, names=NULL, file.out=filepath, open="w")
    
  }
} #0=AA, 1=RNA, 2=DNA, then seqlength

run(1, 600)



