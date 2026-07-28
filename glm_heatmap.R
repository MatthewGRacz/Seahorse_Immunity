pipeline_path <- "Pipeline/AB_62426/"

microbe_attribute_data <- read.csv(paste0(pipeline_path, "GLMOTUSTv2.csv"))
microbe_attribute_data <- column_to_rownames(microbe_attribute_data, "FISH")
microbe_data <- microbe_attribute_data[, grep(x = colnames(microbe_attribute_data), pattern = "^M", value=TRUE)]
#microbe data is the same for the G and H (Microbe vs Supertype and Microbe vs Unique_Recombinant)

get_glm_and_heatmap(pipeline_path, 
                    microbe_attribute_data, 
                    microbe_data, 
                    "S", 
                    "Supertype", 
                    NULL, 
                    NULL, 
                    100, 
                    32, 
                    "JLA Microbe–Supertype Associations for Microbes Present 99+ Times", 
                    "JLA_Supertype_Microbe_WaldZ_Heatmap")

attribute_data <- +(table(analyzed_supertypes_df) > 0) 
attribute_data <- attribute_data[rownames(attribute_data) %in% rownames(microbe_data), ]
colnames(attribute_data) <- sprintf("S%02d", as.numeric(colnames(attribute_data)))
#logical vectors of if a supertype is present or not per individual, just with my supertypes
#only for individuals whose microbes were analyzed
#renames their supertype numbers to the S02d format (eg. S01, S23)

attribute_freqs <- table(sprintf("S%02d", as.numeric(final_dapc$grp)))
#the frequency of each supertype from my analysis
#also renames final_dapc's supertype numbers to the S02d format (eg. S01, S23)

#MGR Microbe vs Supertype Analysis

get_glm_and_heatmap(pipeline_path, 
                    microbe_attribute_data, 
                    microbe_data, 
                    "S", 
                    "Supertype", 
                    attribute_data, 
                    attribute_freqs, 
                    100, 
                    32, 
                    "MGR Microbe–Supertype Associations for Microbes Present 99+ Times", 
                    "MGR_Supertype_Microbe_WaldZ_Heatmap")


microbe_attribute_data <- read.csv(paste0(pipeline_path, "GLMOTUMH-021521.csv"))
microbe_attribute_data <- column_to_rownames(microbe_attribute_data, "FISH")
microbe_data <- microbe_attribute_data[, grep(x = colnames(microbe_attribute_data), pattern = "^M", value=TRUE)]
#microbe data is the same for the G and H (Microbe vs Supertype and Microbe vs Unique_Recombinant)

get_glm_and_heatmap(pipeline_path, 
                    microbe_attribute_data, 
                    microbe_data, 
                    "AB", 
                    "Unique Recombinant", 
                    NULL, 
                    NULL, 
                    0, 
                    32, 
                    "JLA Microbe–Unique Recombinant Associations for All Microbes", 
                    "JLA_UniqueRecombinant_Microbe_WaldZ_Heatmap")


attribute_data <- +(table(recomb_freqs[, c("INDIVIDUAL", "UNIQUE_RECOMB")]) > 0) 
#logical vectors of if a unique_recomb is present or not per individual
#may differ based on what RECCO filters out, which is a bit stochastic, since some recombinants may not be in the
#final kept set of recombinants; the value is similar to JLA's but can differ
#this particular one is already defined for the individuals and their recombinants in JLA's microbe analysis
#also, recomb_freqs already named the unique_recombs with the AB03d format, so no need for renaming

attribute_freqs <- table(recomb_freqs$UNIQUE_RECOMB)
#frequency of a unique_recomb throughout all individuals sampled for microbes

get_glm_and_heatmap(pipeline_path, 
                    microbe_attribute_data, 
                    microbe_data, 
                    "AB", 
                    "Unique Recombinant", 
                    attribute_data, 
                    attribute_freqs, 
                    0, 
                    32, 
                    "MGR Microbe–Unique Recombinant Associations for All Microbes", 
                    "MGR_UniqueRecombinant_Microbe_WaldZ_Heatmap")


get_glm_and_heatmap <- function(pipeline_path, 
                            microbe_attribute_data, 
                            microbe_data,
                            attribute_symbol, 
                            attribute_name, 
                            attribute_data, 
                            attribute_freqs, 
                            min_reads, 
                            num_microbes, 
                            heatmap_title, 
                            heatmap_file_name){
  
  attribute_symbol <- paste0("^", attribute_symbol)
  #for regular expression functions
  
  if(is.null(attribute_data)){
    
    #set as NULL for JLA, since her data contains the logical vectors of presence/absence 
    #of individuals' recombinants in supertypes or unique_recombs or whatever other attribute
    #we can only know the abundance of my attribute, since JLA's frequencies are not available
    #and her logical vector can make it so an individual with 4 recombs all of one attribute 
    #get condensed to a 1 for that attribute, thereby losing data
  
    attribute_data <- microbe_attribute_data[, grep(x = colnames(microbe_attribute_data), pattern = attribute_symbol, value=TRUE)]
    
    #this gets the columns of logical vectors for presence/absence for that attribute by individual
    
  }
  
  used_microbes <- microbe_data[, colSums(microbe_data) >= min_reads, drop = FALSE]
  #filters out microbes with fewer than min_reads number of reads
  
  microbe_attribute_combos <- setNames(
                expand.grid(colnames(used_microbes), 
                colnames(attribute_data), 
                stringsAsFactors = FALSE),
                
    c("Microbe", attribute_name))
  
  
  GLMresults <- apply(microbe_attribute_combos, MARGIN = 1, FUN = function(microbe_attribute) {
    
    attribute <- unname(microbe_attribute[[attribute_name]])
    microbe   <- unname(microbe_attribute["Microbe"])
    
    glm_analysis <- suppressWarnings(coef(summary(glm(attribute_data[, attribute] ~ used_microbes[, microbe], 
                                                      family = quasibinomial(link = "logit")))))
    
    #the actual GLM analysis for each association
    
    data.frame(
      RECOMB_MICROBE = paste0(attribute, microbe), 
      MICROBE = microbe, 
      ATTRIBUTE = attribute, 
      SLOPE = glm_analysis[2, 1],
      SE = glm_analysis[2, 2],
      WALDZ = glm_analysis[2, 1] / glm_analysis[2, 2],
      P = glm_analysis[2, 4],
      stringsAsFactors = FALSE
    )
    
  })
  
  GLMresults <- do.call(rbind, GLMresults)
  
  GLMresults$SEQ_BONFERRONI_P <- p.adjust(GLMresults$P, method = "holm")
  #Sequential Bonferroni correction (Holm-Bonferroni) along the entire dataframe
  
  GLMresults$SIGNIFICANCE <- as.character(symnum(GLMresults$SEQ_BONFERRONI_P, 
                                                   corr = FALSE, 
                                                   na = FALSE, 
                                                   cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                                   symbols = c("***", "**", "*", "")))
  
  GLMresults <- GLMresults[GLMresults$MICROBE %in% names(sort(colSums(used_microbes, na.rm = TRUE), decreasing = TRUE))[1:num_microbes], ]
  
  min_z <- min(GLMresults$WALDZ, na.rm = TRUE)
  max_z <- max(GLMresults$WALDZ, na.rm = TRUE)
  mid_z <- (min_z + max_z)/2
  
  heatmap <- ggplot(GLMresults, aes(x = paste0(ATTRIBUTE, "\n",  attribute_freqs[ATTRIBUTE]), y = MICROBE, fill = WALDZ)) +
    geom_tile(color = "black") +
    geom_text(aes(label = SIGNIFICANCE, color = WALDZ > mid_z), size = 5, vjust = 0.7)  +
    scale_color_manual(values = c("TRUE" = "white", "FALSE" = "black"), guide = "none") +
    scale_fill_gradient2(
      low = "blue", mid="hotpink", high = "#ff0026", 
      midpoint=mid_z,
      limits = c(min_z, max_z),
      breaks = c(min_z, max_z),
      labels = ceiling(c(min_z, max_z)), 
      name = "Wald's Z"
    ) +
    guides(fill = guide_colorbar(frame.color = "black", frame.linewidth = 0.2)) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 8)
    ) +
    labs(title = heatmap_title, x = attribute_name , y = "Microbe")
  
  suppressWarnings(ggsave(paste0(pipeline_path, heatmap_file_name, ".pdf"), 
                          plot = heatmap, 
                          width = 10, 
                          height = 8))
  
}
  
  
  
  
