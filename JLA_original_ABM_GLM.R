#R Script to process MHC/Microbe statistical correlations - ABW (011519) edits JLA (020519) edits (011721)

#library(xlsx) requires rJava package to run that is notoriously problematic on MacOS
#I could not find a fix to get rJava running, so I opted for a Java-free package that writes to xlsx files.

#ABW-011721
#Edited Script to process alleles (AB) instead of supertypes(S)


library(writexl)
library(vegan)

#MHC Allele (IIa/IIb) vs. Microbe abundance (R vs. M)

G <- read.csv(file = "/Users/mattracz/Projects/Wilson_Lab/MatthewRacz-052725/02-Host-Microbe Project/supporting analysis docs/GLM/GLM021721 - Script and Data/GLMOTUMH_021521.csv", header = TRUE)
str(G)
modelAB001M0001 <-glm(formula = G$AB001 ~ G$M0001, family = quasibinomial(link="logit"))
summary(modelAB001M0001)

#MHC Supertype (IIa/IIb) vs. Microbe abundance (S vs. M)

#H <- read.csv(file = "GLMOTUSTv2.csv", header = TRUE)
#str(H)
#modelS00M0001 <-glm(formula = H$S00 ~ H$M0001, family = quasibinomial(link="logit"))
#summary(modelS00M0001)

#MHC Allele Analysis

#Count number of Microbe columns (^M)
NMicrobes<-length(grep(x = colnames(G), pattern = "^M"))

#Count number of Supertype columns (^AB)
NAlleles<-length(grep(x = colnames(G), pattern = "^AB"))

#Total OTUS per fish - 8158 [Data are rarified to 8158]

FishOTUS<-apply(G[2:(NMicrobes+1)], 1, sum)


#Extract microbial OTUS present >99X, retain FISH ID
Gsubset<-G[c(rep(TRUE, 1), colSums(G[2:(NMicrobes+1)]) > 99)]

#Extract Microbe names
Microbes<-grep(x = names(Gsubset), pattern = "^M",value=TRUE)

#Extract Allele names
Alleles<-grep(x = names(Gsubset), pattern = "^AB",value=TRUE)

#Create a dataframe to store results
AlleleMicrobeAnalysis<-data.frame()

#Microbes <- c("M0104", "M0102", "M0064", "M0044", "M0037", "M0030", "M0026", "M0025", "M0024", "M0023", "M0022", "M0021", "M0020", "M0019", "M0018", "M0017", "M0016", "M0015", "M0014", "M0013", "M0012", "M0011", "M0010", "M0009", "M0008", "M0007", "M0006", "M0005", "M0004", "M0003", "M0002", "M0001")

#GLM loop - store intercept, slope, se, p.value
for (valM in Microbes) {
  for (valAB in Alleles) {
    M<-valM
    AB<-valAB
    ABM<-paste(valAB,valM,sep="")
    
    fit <-glm(formula = Gsubset[,valAB] ~ Gsubset[,valM], family = quasibinomial(link="logit"))
    
    ## capture summary stats
    intercept <- coef(summary(fit))[1]
    slope <- coef(summary(fit))[2]
    se <- coef(summary(fit))[4]
    WaldZ <- slope/se
    p.value <- coef(summary(fit))[8]
    
    # create temporary data frame
    df <- data.frame(AlleleMicrobe = ABM, intercept = intercept,
                     slope = slope, se = se, WALDZ = WaldZ, p.value = p.value, MICROBE = valM, ATTRIBUTE = valAB, stringsAsFactors = F)
    
    
    # bind rows of temporary data frame to the results data frame
    AlleleMicrobeAnalysis <- rbind(AlleleMicrobeAnalysis, df)
  }
}     

AlleleMicrobeAnalysis <- AlleleMicrobeAnalysis[AlleleMicrobeAnalysis$MICROBE %in% names(sort(colSums(Gsubset[, Microbes], na.rm = TRUE), decreasing = TRUE))[1:32], ]


#Write data to xlsx
#write.xlsx(AlleleMicrobeAnalysis, "AlleleMicrobeAnalysis.xlsx") 
#using writexl

#Relative microbe data

Grelative<-G[c(rep(TRUE, 1), colSums(G[2:(NMicrobes+1)]) > -1)]
rownames(Grelative) <- Grelative[,1]
Grelative[,1] <- NULL

#Total microbe abundances

apply(Grelative[1:length(Microbes)], 2,sum)

#Calculate microbe relative abundances

Grelative <- cbind(decostand(Grelative[1:length(Microbes)], method = "total"), Grelative[(length(Microbes)+1):length(Grelative)])


#Create a dataframe to store results
AlleleRelativeMicrobeAnalysis<-data.frame()

#GLM loop - store intercept, slope, p.value
for (valM in Microbes) {
  for (valAB in Alleles) {
    M<-valM
    AB<-valAB
    ABM<-paste(valAB,valM,sep="")
    
    fit <-glm(formula = Grelative[,valAB] ~ Grelative[,valM], family = quasibinomial(link="logit"))
    
    ## capture summary stats
    intercept <- coef(summary(fit))[1]
    slope <- coef(summary(fit))[2]
    se <- coef(summary(fit))[4]
    WaldZ <- slope/se
    p.value <- coef(summary(fit))[8]
    
    # create temporary data frame
    df <- data.frame(AlleleMicrobe = ABM, intercept = intercept,
                     slope = slope, WALDZ = WaldZ, p.value = p.value, MICROBE = valM, ATTRIBUTE = valAB, stringsAsFactors = F)
    
    
    # bind rows of temporary data frame to the results data frame
    AlleleRelativeMicrobeAnalysis <- rbind(AlleleRelativeMicrobeAnalysis, df)
  }
}

AlleleRelativeMicrobeAnalysis <- AlleleRelativeMicrobeAnalysis[AlleleRelativeMicrobeAnalysis$MICROBE %in% names(sort(colSums(Grelative[, Microbes], na.rm = TRUE), decreasing = TRUE))[1:32], ]

#Write data to xlsx
#write.xlsx(AlleleMicrobeAnalysis, "AlleleMicrobeAnalysis.xlsx")
#using writexl

AlleleMicrobeAnalysis$p.value <- p.adjust(AlleleMicrobeAnalysis$p.value, method = "holm")

AlleleMicrobeAnalysis$SIGNIFICANCE <- as.character(symnum(AlleleMicrobeAnalysis$p.value, 
                                                          corr = FALSE, 
                                                          na = FALSE, 
                                                          cutpoints = c(0, 0.00001, 1), 
                                                          symbols = c("**", "")))

write.csv(AlleleMicrobeAnalysis, file="/Users/mattracz/Projects/Wilson_Lab/Pipeline/CSV.csv")

#get the top X most abundant microbes from the GLM results, to show for the heatmap

min_z <- min(AlleleMicrobeAnalysis$WALDZ, na.rm = TRUE)
max_z <- max(AlleleMicrobeAnalysis$WALDZ, na.rm = TRUE)
mid_z <- 0

abs_heatmap <- ggplot(AlleleMicrobeAnalysis, aes(x = ATTRIBUTE, y = MICROBE, fill = WALDZ)) +
  geom_tile(color = "black", size=0.4) +
  geom_text(aes(label = SIGNIFICANCE, color = WALDZ > mid_z), size = 5, vjust = 0.7)  +
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
    axis.text.x = element_text(size = 4), 
    axis.text.y = element_text(size = 10), 
    plot.title = element_text(size = 20 * 0.5, face = "bold"), 
    axis.title.x = element_text(size = 0.5 * 18, face = "bold"),
    axis.title.y = element_text(size = 0.5 * 50, face = "bold")
  ) +
  labs(title = "JLA original Absolute AB-Microbe Heatmap With Bonferroni Correction", x = "AB" , y = "Microbe")

AlleleRelativeMicrobeAnalysis$p.value <- p.adjust(AlleleRelativeMicrobeAnalysis$p.value, method = "holm")

AlleleRelativeMicrobeAnalysis$SIGNIFICANCE <- as.character(symnum(AlleleRelativeMicrobeAnalysis$p.value, 
                                                                  corr = FALSE, 
                                                                  na = FALSE, 
                                                                  cutpoints = c(0, 0.00001, 1), 
                                                                  symbols = c("**", "")))

#get the top X most abundant microbes from the GLM results, to show for the heatmap

min_z <- min(AlleleRelativeMicrobeAnalysis$WALDZ, na.rm = TRUE)
max_z <- max(AlleleRelativeMicrobeAnalysis$WALDZ, na.rm = TRUE)
mid_z <- 0

rel_heatmap <- ggplot(AlleleRelativeMicrobeAnalysis, aes(x = ATTRIBUTE, y = MICROBE, fill = WALDZ)) +
  geom_tile(color = "black", size=0.4) +
  geom_text(aes(label = SIGNIFICANCE, color = WALDZ > mid_z), size = 5, vjust = 0.7)  +
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
    axis.text.x = element_text(size = 4), 
    axis.text.y = element_text(size = 10), 
    plot.title = element_text(size = 20 * 0.5, face = "bold"), 
    axis.title.x = element_text(size = 0.5 * 18, face = "bold"),
    axis.title.y = element_text(size = 0.5 * 50, face = "bold")
  ) +
  labs(title = "JLA original Relative AB-Microbe Heatmap With Bonferroni Correction", x = "AB" , y = "Microbe")

dev.new()
print(abs_heatmap)
print(rel_heatmap)

#suppressWarnings(ggsave(filename = "/Users/mattracz/Projects/Wilson_Lab/Pipeline/JLA_original_Abs_AB_Microbes_WithoutBC", plot = abs_heatmap, limitsize = FALSE))
#suppressWarnings(ggsave(filename = "/Users/mattracz/Projects/Wilson_Lab/Pipeline/JLA_original_Rel_AB_Microbes_WithoutBC", plot = rel_heatmap, limitsize = FALSE))
