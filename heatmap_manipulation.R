
#pipeline_path <- "Pipeline/AB_NOCLONES/"
pipeline_path <- "Pipeline/AB_62426/"
#pipeline_path <- "Pipeline/AB_72726/"


zscores <- read.csv(paste0(pipeline_path, "zscores.csv"))
zscores <<- zscores

JLA_microbes <- c('M0008', 'M0003', 'M0006', 'M0016', 'M0044', 'M0018', 'M0061', 'M0081', 'M0031', 'M0014', 'M0012', 'M0007', 'M0015', 'M0019', 'M0011', 'M0005', 'M0010', 'M0030', 'M0041', 'M0085', 'M0029', 'M0020', 'M0048', 'M0013', 'M0002', 'M0004', 'M0023', 'M0032', 'M0064', 'M0001', 'M0017', 'M0009')
#microbes JLA used for Ch 3 heatmap, in the order she has them from bottom up (how the heatmap automatically names the y-axis)


GS_OG_jumper_df <- read.csv("Pipeline/AB_NOCLONES/jumpers_dataframe.csv")

GS_OG_recombs_data <- GS_OG_jumper_df[GS_OG_jumper_df$INDIVIDUAL %in% ordered_names, ]

GS_OG_individuals_data <- table(GS_OG_recombs_data$INDIVIDUAL,factor(GS_OG_recombs_data$BEST, levels = 1:17))

colnames(GS_OG_individuals_data) <- sprintf("S%02d", 0:16)
#uses the supertype names JLA had in her presence/absence matrix

GS_OG_individuals_data <- GS_OG_individuals_data[ordered_names,]
GS_OG_individuals_data <- as.matrix(GS_OG_individuals_data)
#set colnames that JLA had and put individuals in her order

#Uses the original GS data that JLA had, normalized via the Hungarian Algorithm. In sum, this can 
#backtrace her supertype assignments using her presence/absence matrix, and then I can use that
#as a translator/dictionary for newer datasets; new data STs --> old data STs --> JLA STs
#this standardizes every heatmap and hence comparable 1:1 with any filterings/modifications to them

GLM_recomb_jumper_df <- read.csv(paste0(pipeline_path, "jumpers_dataframe.csv"))

GLM_recombs_data <- GLM_recomb_jumper_df[GLM_recomb_jumper_df$INDIVIDUAL %in% ordered_names, ]

GLM_individuals_data <- table(GLM_recombs_data$INDIVIDUAL,factor(GLM_recombs_data$BEST, levels = 1:17))

colnames(GLM_individuals_data) <- sprintf("S%02d", 0:16)

GLM_individuals_data <- GLM_individuals_data[ordered_names,]
GLM_individuals_data <- as.matrix(GLM_individuals_data)


microbe_attribute_data <- read.csv(paste0("/Users/mattracz/Projects/Wilson_Lab/Pipeline/GLMOTUSTv2.csv"))
microbe_attribute_data <- column_to_rownames(microbe_attribute_data, "FISH")
JLA_matrix <- microbe_attribute_data[!grepl("^M", colnames(microbe_attribute_data))]
microbe_attribute_data <- microbe_attribute_data[!grepl("^S", colnames(microbe_attribute_data))]

JLA_matrix[,setdiff(colnames(GLM_individuals_data), colnames(JLA_matrix))] <- 0

JLA_matrix <- JLA_matrix[, colnames(GLM_individuals_data)]
JLA_matrix <- as.matrix(JLA_matrix)
#pad nonexistent supertypes with 0s and arranges supertypes in the correct order

OG_JLA_overlap_matrix <- t(GS_OG_individuals_data) %*% JLA_matrix
#dot product of JLA's presence/absence matrix and my OG supertype assignments
#this builds a dictionary where I can first translate any supertype assignments to the
#ones used in the GS that JLA originally had, then I can take those assignments and translate
#them to JLA's used supertypes, hence standardizing the supertypes of all heatmaps used in the future

OG_to_JLA_dictionary <- setNames(colnames(JLA_matrix)[as.numeric(solve_LSAP(OG_JLA_overlap_matrix, maximum = TRUE))], colnames(GS_OG_individuals_data))


OG_GLM_overlap_matrix <- t(GLM_individuals_data) %*% GS_OG_individuals_data

GLM_to_OG_dictionary <- setNames(colnames(GS_OG_individuals_data)[as.numeric(solve_LSAP(OG_GLM_overlap_matrix, maximum = TRUE))],colnames(GLM_individuals_data))


GLM_to_JLA_dictionary <- setNames(OG_to_JLA_dictionary[GLM_to_OG_dictionary], names(GLM_to_OG_dictionary))



colnames(GLM_individuals_data) <- GLM_to_JLA_dictionary[colnames(GLM_individuals_data)] 
#align column names as they're fed into the dictionary and translated

GLM_individuals_data <- GLM_individuals_data[, colnames(JLA_matrix)]

GLM_individuals_data <- +(GLM_individuals_data[, colnames(JLA_matrix)] > 0)
#presence/absence vector of supertype presence per recombinants of those individuals; 
#the GLM is quasibinomial and hence operates based off these, instead of by quantity 

GLM_individuals_data <- unclass(GLM_individuals_data)
class(GLM_individuals_data) <- "matrix"

#casts it as a matrix, rids of all table features, to support different GLM types and
#prevent errors when cbinding onto the microbe_attribute_data

microbe_attribute_data <- cbind(microbe_attribute_data, GLM_individuals_data[rownames(microbe_attribute_data), ])
#add the translated presence/absence vector to the microbe_attribute data, which is then fed into the GLM function



get_glm_and_heatmap <- function(pipeline_path, 
                                microbe_attribute_data, 
                                JLA_microbes,
                                attribute_symbol, 
                                attribute_name, 
                                attribute_data, 
                                attribute_freqs, 
                                heatmap_title, 
                                heatmap_file_name, 
                                x_scale,
                                y_scale){
  
  attribute_symbol <- paste0("^", attribute_symbol)
  #for regular expression (regex) functions
  
  microbe_data <- microbe_attribute_data[, colnames(microbe_attribute_data) %in% JLA_microbes]
  #microbe data is the same for the G and H (Microbe vs Supertype and Microbe vs Unique_Recombinant)
  
  attribute_data <- attribute_data[rownames(microbe_data), , drop = FALSE]
  #puts it in the same order as JLA
  
  if(is.null(attribute_data)){
    
    #set as NULL for JLA, since her data contains the logical vectors of presence/absence 
    #of individuals' recombinants in supertypes or unique_recombs or whatever other attribute
    #we can only know the abundance of my attribute, since JLA's frequencies are not available
    #and her logical vector can make it so an individual with 4 recombs all of one attribute 
    #get condensed to a 1 for that attribute, thereby losing data
    
    attribute_data <- microbe_attribute_data[, grep(x = colnames(microbe_attribute_data), pattern = attribute_symbol, value=TRUE)]
    
    #this gets the columns of logical vectors for presence/absence for that attribute by individual
    
  }
  
  microbe_attribute_combos <- setNames(
    expand.grid(colnames(microbe_data), 
                colnames(attribute_data), 
                stringsAsFactors = FALSE),
    
    c("Microbe", attribute_name))
  
  #creates every Microbe-Attribute combo to run a GLM on
  
  cat("\nCalcuating GLM... This may take a moment...\n")
  
  cl <- makeCluster(detectCores() - 1)
  #uses all except 1 computer CPU core to run, since can take up to 8 hours!
  
  clusterExport(cl, c("attribute_data", "microbe_data", "attribute_name"), envir = environment())
  #exports relevant data to each subordinate R session that runs in parallel (parApply next line),
  #each of which creates its own mini-environment
  
  GLMresults <- suppressWarnings(parApply(cl, microbe_attribute_combos, MARGIN = 1, FUN = function(microbe_attribute) {
    
    attribute <- unname(microbe_attribute[[attribute_name]])
    microbe   <- unname(microbe_attribute["Microbe"])
    
    #glm_analysis <- suppressWarnings(coef(summary(glm(attribute_data[, attribute] ~ microbe_data[, microbe], family = quasipoisson(link = "log")))))
    
    glm_analysis <- suppressWarnings(coef(summary(glm(attribute_data[, attribute] ~ microbe_data[, microbe], family = quasibinomial(link = "logit")))))
    
    
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
    
  }))
  
  stopCluster(cl)
  #return to normal CPU allocation and function after the parallel computations are done
  
  GLMresults <- do.call(rbind, GLMresults)
  
  GLMresults$SEQ_BONFERRONI_P <- p.adjust(GLMresults$P, method = "holm")
  #Sequential Bonferroni correction (Holm-Bonferroni) along the entire dataframe
  
  GLMresults$SIGNIFICANCE <- as.character(symnum(GLMresults$SEQ_BONFERRONI_P, 
                                                 corr = FALSE, 
                                                 na = FALSE, 
                                                 cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                                 symbols = c("***", "**", "*", "")))
  
  #the corrected p value per GLM analysis gets marked with stars if p <= certain thresholds
  
  GLMresults$MICROBE <- factor(GLMresults$MICROBE, levels = JLA_microbes)
  #keeps the microbes in the order they're in
  
  cat("\nGLM calculations complete!\n")
  
  min_z <- min(GLMresults$WALDZ, na.rm = TRUE)
  max_z <- max(GLMresults$WALDZ, na.rm = TRUE)
  mid_z <- (min_z + max_z)/2
  
  heatmap <- suppressWarnings(ggplot(GLMresults, aes(x = paste0(ATTRIBUTE, "\n",  attribute_freqs[ATTRIBUTE]), y = MICROBE, fill = WALDZ)) +
                                geom_tile(color = "black", size=0.4) +
                                geom_text(aes(label = SIGNIFICANCE, color = WALDZ > mid_z), size = 4, vjust = 0.7)  +
                                scale_color_manual(values = c("TRUE" = "white", "FALSE" = "black"), guide = "none") +
                                scale_fill_gradient2(
                                  low = "white", mid="#A6A6A6", high = "black", 
                                  midpoint=mid_z,
                                  limits = c(min_z, max_z),
                                  breaks = c(min_z, max_z),
                                  labels = ceiling(c(min_z, max_z)), 
                                  name = "Wald's Z"
                                ) +
                                guides(fill = guide_colorbar(frame.color = "black", frame.linewidth = 0.2)) +
                                theme_minimal() +
                                scale_x_discrete(position = "top") +
                                theme(
                                  axis.text.x = element_text(size = x_scale * 13), 
                                  axis.text.y = element_text(size = y_scale * 30), 
                                  plot.title = element_text(size = 70 * x_scale, face = "bold"), 
                                  axis.title.x = element_text(size = x_scale * 50, face = "bold"),
                                  axis.title.y = element_text(size = y_scale * 50, face = "bold")
                                ) +
                                labs(title = heatmap_title, x = attribute_name , y = "Microbe"))
  
  #suppressWarnings(ggsave(paste0(pipeline_path, heatmap_file_name, ".pdf"), plot = heatmap, width = x_scale * length(unique(GLMresults$ATTRIBUTE)), height = num_microbes * y_scale, limitsize = FALSE))
  
  print(heatmap)
  
  return(GLMresults)
  
}



JLA_MST_data <- get_glm_and_heatmap(pipeline_path, 
                                    microbe_attribute_data, 
                                    JLA_microbes,
                                    "S", 
                                    "Supertype", 
                                    NULL, 
                                    NULL, 
                                    "JLA Microbe–Supertype Associations for Microbes Present 99+ Times", 
                                    "JLA_Supertype_Microbe_WaldZ_Heatmap",
                                    0.3,
                                    0.3)

#write_csv(JLA_MST_data, paste0(pipeline_path, "Microbe_Supertype_GLM_data_QUASIPOISSON.csv"))
write_csv(JLA_MST_data, paste0(pipeline_path, "Microbe_Supertype_GLM_data.csv"))





