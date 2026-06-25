origDNA <- "ATCCACACAGAAATGCAAATGATTGGCTGCTCTGAGACGGAARGAGCATATGACRTCTCTCTGGACGAGGAGGACGTGATCGTTGCCGATTTCATCGCAGAGRGGGTTGTGGATAMACAGCCCAGCTTTGTGGATCACGTCAGCTTTGSCAAAGATRYTTTTCAGGTGGCAATGGGGAATCTAATCACATGCAGAAKAAACCTCAACATAACCCAGAAAGCATTGAAGGATCAGCCACTAGAATCTATGGCTACCTGTKKCATRCGGACKKTGGTTGTGTGTTCAACTCGAGTGACCYGAATGACATCGAGTACTTCSAGATTTACAACTACAACAAACTGAAGCTTTTCCGCTTCAGCAGCACTTTGGATAAGTACGTCGGCTACACMGAGTTTGGCATCAAGCAGGCTACCGCCTTCAACAACAACAAAGACYWCATCGCCGACGTCAGAGCCAGGAAAGAAYACATTTGTTTAMACAATATTAAGAWTGACTACRAAAGTGCGCTCACCAAGTCAG"

origDNA <- toupper(origDNA)
origDNA <- strsplit(origDNA, "")[[1]]

first_allele <- "ATCCACACAGAAATGCAAATGATTGGCTGCTCTGAGACGGAAAGAGCATATGACATCTCTCTGGACGAGGAGGACGTGATCGTTGCCGATTTCATCGCAGAGAGGGTTGTGGATACACAGCCCAGCTTTGTGGATCACGTCAGCTTTGGCAAAGATATTTTTCAGGTGGCAATGGGGAATCTAATCACATGCAGAAGAAACCTCAACATAACCCAGAAAGCATTGAAGGATCAGCCACTAGAATCTATGGCTACCTGTTTCATACGGACTTTGGTTGTGTGTTCAACTCGAGTGACCCGAATGACATCGAGTACTTCGAGATTTACAACTACAACAAACTGAAGCTTTTCCGCTTCAGCAGCACTTTGGATAAGTACGTCGGCTACACCGAGTTTGGCATCAAGCAGGCTACCGCCTTCAACAACAACAAAGACTTCATCGCCGACGTCAGAGCCAGGAAAGAATACATTTGTTTAAACAATATTAAGATTGACTACGAAAGTGCGCTCACCAAGTCAG"

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


