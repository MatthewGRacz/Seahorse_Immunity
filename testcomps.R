


unique_individuals <- unique(allele_seq_df$INDIVIDUAL)

compare_alleles <- function(indv, ab_col, a_b_col, AB_sequences, A_B_sequences, init_index){
  
  indv_data <- allele_seq_df[allele_seq_df$INDIVIDUAL == indv, ]
  
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

AB_sequences <- suppressWarnings(toupper(read.fasta('/Users/mattracz/Projects/Wilson_Lab/PHASED/p=0.5/NO_CLONES/PHASE_AB_NOCLONES/unphased.fasta', as.string = TRUE)))

sapply(unique_individuals, 
       compare_alleles, 
       ab_col = "AB_ALPHA", 
       a_b_col = "ALPHA", 
       AB_sequences = allele_seq_df[["UNPHASED_AB"]], 
       A_B_sequences = allele_seq_df[["UNPHASED_ALPHA"]], 
       init_index=0)


cat("\n\nBeta comparisons:\n\n")

sapply(unique_individuals, 
       compare_alleles, 
       ab_col = "AB_BETA", 
       a_b_col = "BETA", 
       AB_sequences = allele_seq_df[["UNPHASED_AB"]], 
       A_B_sequences = allele_seq_df[["UNPHASED_BETA"]], 
       init_index=246)


cat("\nhello")
