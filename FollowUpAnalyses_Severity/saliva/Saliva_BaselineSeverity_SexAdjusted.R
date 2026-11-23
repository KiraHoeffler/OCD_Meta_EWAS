
################################################################################
# SETUP
################################################################################

# LOAD LIBRARIES
library(DBI, lib.loc = "Z:/Bioconductorpackages/319")
library(tidyselect, lib.loc = "Z:/Bioconductorpackages/319")
library(shiny, lib.loc = "Z:/Bioconductorpackages/319")
library(readxl)
library(qqman)
library(ggplot2)
library(Matrix)
library(tictoc)
library(writexl)
library(dplyr)
library(tidyverse)
library(sva)
library(car)
library(limma)
library(ENmix)
library(patchwork)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(Gviz)
library(magick)
library(microbenchmark)
library(COMBAT)
library(pscl)

# SET PATHS
dir_gen <- "S://Project/WP-epigenetics/11_BaselineSeverity/" #working directory
beta_path <- "S://Project/WP-epigenetics/04_Pipeline_combinedPhases/output/RData/beta_final_autosomes.RData"
Mvalues_path <- "S://Project/WP-epigenetics/04_Pipeline_combinedPhases/output/RData/Mvalues_final_autosomes.RData"
anno_path <- "S://Project/WP-epigenetics/02_Import/Updated_annotation/UPDATED_ANNOTATION.RData"
ctrl_path <- "S://Project/WP-epigenetics/04_Pipeline_combinedPhases/output/RData/Positive_ctrlprobe_intensities.RData"
SNP_path <- "S://Project/WP-epigenetics/04_Pipeline_combinedPhases/output/RData/comb_SNPs.RData"
samplesheet_path <- "S://Project/WP-epigenetics/04_Pipeline_combinedPhases/samplesheet_clinic.xlsx"

################################################################################
# IMPORT
################################################################################

# SET WORKING DIRECTORY
setwd(dir_gen)

# LOAD SAMPLESHEET & EUROPEAN INFO
samplesheet <- as.data.frame(read_excel(samplesheet_path))

# LOAD M VALUES
load(Mvalues_path)
Mvalues <- Mvalues[row.names(Mvalues) != "cg05575921", ] # remove cpg that is used to adjust for smoking

# CTRL PROBES
load(ctrl_path) 

# SNP PROBES
load(SNP_path) 

# BETA VALUES
load(beta_path)

# anno_updTATION
load(anno_path)

# GGPLOT2 THEME
load("output_CrossPhase/RData/theme.RData")
load("output_CrossPhase/RData/theme_transparent.Rdata")

# IMPORT FUNCTIONS
source(file = "../05_Response/0_General/Functions/DMR_functions.R")
source(file = "../05_Response/0_General/Functions/PCAs_for_model.R")
source(file = "../05_Response/0_General/Functions/lambda.R")

# MAKE FOLDERS
if (!dir.exists("output_CrossPhase/Figures")) dir.create("output_CrossPhase/Figures")
if (!dir.exists("output_CrossPhase/Figures/DMRs")) dir.create("output_CrossPhase/Figures/DMRs")
if (!dir.exists("output_CrossPhase/Figures/DMRs/raw")) dir.create("output_CrossPhase/Figures/DMRs/raw")
if (!dir.exists("output_CrossPhase/Figures/DMRs/combined")) dir.create("output_CrossPhase/Figures/DMRs/combined")
if (!dir.exists("output_CrossPhase/Figures/DMRs/Comorbidities")) dir.create("output_CrossPhase/Figures/DMRs/Comorbidities")
if (!dir.exists("output_CrossPhase/Figures/DMRs/Medicine")) dir.create("output_CrossPhase/Figures/DMRs/Medicine")
if (!dir.exists("output_CrossPhase/Figures/DMRs/Controls")) dir.create("output_CrossPhase/Figures/DMRs/Controls")
if (!dir.exists("output_CrossPhase/Figures/DMRs/mQTL")) dir.create("output_CrossPhase/Figures/DMRs/mQTL")
if (!dir.exists("output_CrossPhase/Figures/DMRs/mQTL/raw")) dir.create("output_CrossPhase/Figures/DMRs/mQTL/raw")
if (!dir.exists("output_CrossPhase/Figures/DMRs/mQTL/combined")) dir.create("output_CrossPhase/Figures/DMRs/mQTL/combined")
if (!dir.exists("output_CrossPhase/Figures/DMRs/mQTL_DNAm")) dir.create("output_CrossPhase/Figures/DMRs/mQTL_DNAm")
if (!dir.exists("output_CrossPhase/Tables")) dir.create("output_CrossPhase/Tables")
if (!dir.exists("output_CrossPhase/RData")) dir.create("output_CrossPhase/RData")

################################################################################
# GENERAL
################################################################################

# EXCLUDE THOSE THAT DO NOT HAVE A severity SCORE AT 3 MONTHS
samples_resp <- samplesheet[which(samplesheet$Time_point == "D1" & !is.na(samplesheet$Pre_YB_Sum)), ]
nrow(samples_resp) #898
write_xlsx(samples_resp, path = "output_CrossPhase/Tables/Samplesheet_severity.xlsx")

# FILTER M VALUES ON DAY1 SAMPLES WITH severity
Mvalues_filt <- Mvalues[, samples_resp$Basename]


################################################################################
# MODEL PREPARATION
################################################################################

# PCs FOR MODEL
ctrlprobe_PCAscores <- ctrl_probe_PCA(ctrl, samples_resp)
Ctrl_PC1 = ctrlprobe_PCAscores$PC1
Ctrl_PC2 = ctrlprobe_PCAscores$PC2
Ctrl_PC3 = ctrlprobe_PCAscores$PC3
Ctrl_PC4 = ctrlprobe_PCAscores$PC4
Ctrl_PC5 = ctrlprobe_PCAscores$PC5
Ctrl_PC6 = ctrlprobe_PCAscores$PC6
Ctrl_PC7 = ctrlprobe_PCAscores$PC7
Ctrl_PC8 = ctrlprobe_PCAscores$PC8
Ctrl_PC9 = ctrlprobe_PCAscores$PC9
Ctrl_PC10 = ctrlprobe_PCAscores$PC10
Ctrl_PC11 = ctrlprobe_PCAscores$PC11
Ctrl_PC12 = ctrlprobe_PCAscores$PC12
Ctrl_PC13 = ctrlprobe_PCAscores$PC13
Ctrl_PC14 = ctrlprobe_PCAscores$PC14
Ctrl_PC15 = ctrlprobe_PCAscores$PC15

cellcounts_PCAscores <- cell_type_PCA(samples_resp, c("Epi", "Fib", "comb_ICs"))
cell_PC1 = cellcounts_PCAscores$PC1
cell_PC2 = cellcounts_PCAscores$PC2

snpsprobe_PCAscores <- ancestry_PCA(comb_SNPs, samples_resp)
snps_PC1 = snpsprobe_PCAscores$PC1
snps_PC2 = snpsprobe_PCAscores$PC2
snps_PC3 = snpsprobe_PCAscores$PC3
snps_PC4 = snpsprobe_PCAscores$PC4
snps_PC5 = snpsprobe_PCAscores$PC5
snps_PC6 = snpsprobe_PCAscores$PC6
snps_PC7 = snpsprobe_PCAscores$PC7
snps_PC8 = snpsprobe_PCAscores$PC8
snps_PC9 = snpsprobe_PCAscores$PC9
snps_PC10 = snpsprobe_PCAscores$PC10


rm(ctrlprobe_PCAscores, cellcounts_PCAscores, snpsprobe_PCAscores)

# VARIABLES FOR MODEL
severity = scale(as.numeric(samples_resp$Pre_YB_Sum))
Age = scale(as.numeric(samples_resp$Age))
smoking = scale(as.numeric(samples_resp$smoking))
Sex <- as.factor(samples_resp$Sex)
Phase <- as.factor(samples_resp$Phase)
Collection_Year <- as.factor(samples_resp$Collection_year)

# Check for multi-colinearity
vif_obj <- as.data.frame(vif(lm(Mvalues_filt[1, ] ~ severity + Sex + Age + smoking + Phase + Collection_Year + 
                                  Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 + Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10 +
                                  Ctrl_PC11 +  Ctrl_PC12 +  Ctrl_PC13 +  Ctrl_PC14 +  Ctrl_PC15 +
                                  cell_PC1 + cell_PC2 +
                                  snps_PC1 + snps_PC2 + snps_PC3 + snps_PC4 + snps_PC5 + snps_PC6 + snps_PC7 + snps_PC8 + snps_PC9 + snps_PC10)))
vif_obj$Factor <- rownames(vif_obj)
colnames(vif_obj) <- c("vif", "factor")
vif_obj
write_xlsx(vif_obj, path = "output_CrossPhase/Tables/VIF.xlsx")


#after excluding ctrl_pc4 because of multi-colinearity
vif_obj <- as.data.frame(vif(lm(Mvalues_filt[1, ] ~ severity + Sex + Age + smoking + Collection_Year + 
                                  Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC5 + Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10 +
                                  Ctrl_PC11 +  Ctrl_PC12 +  Ctrl_PC13 +  Ctrl_PC14 +  Ctrl_PC15 +
                                  cell_PC1 + cell_PC2 +
                                  snps_PC1 + snps_PC2 + snps_PC3 + snps_PC4 + snps_PC5 + snps_PC6 + snps_PC7 + snps_PC8 + snps_PC9 + snps_PC10)))
vif_obj$Factor <- rownames(vif_obj)
colnames(vif_obj) <- c("vif", "factor")
vif_obj
write_xlsx(vif_obj, path = "output_CrossPhase/Tables/VIF2.xlsx")


################################################################################
# RUN DIFFERENT MODELS
################################################################################

# TABLE TO SAVE ALL LAMBDA VALUES
lambda_table <- data.frame(Model = c("2Ctrl_2Anc", "5Ctrl_2Anc", "10Ctrl_2Anc", "15Ctrl_2Anc",
                                     "2Ctrl_5Anc", "5Ctrl_5Anc", "10Ctrl_5Anc", "15Ctrl_5Anc",
                                     "2Ctrl_10Anc", "5Ctrl_10Anc", "10Ctrl_10Anc", "15Ctrl_10Anc",
                                     "2Ctrl_1Anc", "5Ctrl_1Anc", "10Ctrl_1Anc", "15Ctrl_1Anc"), lambda = NA)


### 2 Ctrl PCs, 2 Anc PCs ###
design <- model.matrix(Mvalues_filt[1, ] ~ severity + Sex + Age + smoking + Collection_Year + cell_PC1 + cell_PC2 + 
                         Ctrl_PC1 + Ctrl_PC2 +
                         snps_PC1 + snps_PC2)
fit <- lmFit(Mvalues_filt, design)
fit <- eBayes(fit)

DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)

pdf(file = "output_CrossPhase/Figures/QQplot_Baselineseverity_2CtrlPCs_2AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output_CrossPhase/Figures/QQplot_Baselineseverity_2CtrlPCs_2AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[1,2] <- lambda(DMPs$P.Value)


### 5 Ctrl PCs, 2 Anc PCs ###
design <- model.matrix(Mvalues_filt[1, ] ~ severity + Sex + Age + smoking + Collection_Year + cell_PC1 + cell_PC2 +
                         Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC5 +
                         snps_PC1 + snps_PC2)
fit <- lmFit(Mvalues_filt, design)
fit <- eBayes(fit)

DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)

pdf(file = "output_CrossPhase/Figures/QQplot_Baselineseverity_5CtrlPCs_2AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output_CrossPhase/Figures/QQplot_Baselineseverity_5CtrlPCs_2AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[2,2] <- lambda(DMPs$P.Value)


### 10 Ctrl PCs, 2 Anc PCs ###
design <- model.matrix(Mvalues_filt[1, ] ~ severity + Sex + Age + smoking + Collection_Year + cell_PC1 + cell_PC2 +
                         Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC5 + Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10 +
                         snps_PC1 + snps_PC2)
fit <- lmFit(Mvalues_filt, design)
fit <- eBayes(fit)

DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)

pdf(file = "output_CrossPhase/Figures/QQplot_Baselineseverity_10CtrlPCs_2AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output_CrossPhase/Figures/QQplot_Baselineseverity_10CtrlPCs_2AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[3,2] <- lambda(DMPs$P.Value)


### 15 Ctrl PCs, 2 Anc PCs ###
design <- model.matrix(Mvalues_filt[1, ] ~ severity + Sex + Age + smoking + Collection_Year + cell_PC1 + cell_PC2 + 
                         Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC5 + Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10 +
                         Ctrl_PC11 + Ctrl_PC12 + Ctrl_PC13 + Ctrl_PC14 + Ctrl_PC15 +
                         snps_PC1 + snps_PC2)
fit <- lmFit(Mvalues_filt, design)
fit <- eBayes(fit)

DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)

pdf(file = "output_CrossPhase/Figures/QQplot_Baselineseverity_15CtrlPCs_2AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output_CrossPhase/Figures/QQplot_Baselineseverity_15CtrlPCs_2AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[4,2] <- lambda(DMPs$P.Value)



### 2 Ctrl PCs, 5 Anc PCs ###
design <- model.matrix(Mvalues_filt[1, ] ~ severity + Sex + Age + smoking + Collection_Year + cell_PC1 + cell_PC2 + 
                         Ctrl_PC1 + Ctrl_PC2 + 
                         snps_PC1 + snps_PC2 + snps_PC3 + snps_PC4 + snps_PC5)
fit <- lmFit(Mvalues_filt, design)
fit <- eBayes(fit)

DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)

pdf(file = "output_CrossPhase/Figures/QQplot_Baselineseverity_2CtrlPCs_5AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output_CrossPhase/Figures/QQplot_Baselineseverity_2CtrlPCs_5AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[5,2] <- lambda(DMPs$P.Value)


### 5 Ctrl PCs, 5 Anc PCs ###
design <- model.matrix(Mvalues_filt[1, ] ~ severity + Sex + Age + smoking + Collection_Year + cell_PC1 + cell_PC2 + 
                         Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC5 +
                         snps_PC1 + snps_PC2 + snps_PC3 + snps_PC4 + snps_PC5)
fit <- lmFit(Mvalues_filt, design)
fit <- eBayes(fit)

DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)

pdf(file = "output_CrossPhase/Figures/QQplot_Baselineseverity_5CtrlPCs_5AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output_CrossPhase/Figures/QQplot_Baselineseverity_5CtrlPCs_5AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[6,2] <- lambda(DMPs$P.Value)


### 10 Ctrl PCs, 5 Anc PCs ###
design <- model.matrix(Mvalues_filt[1, ] ~ severity + Sex + Age + smoking + Collection_Year + cell_PC1 + cell_PC2 + 
                         Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC5 + Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10 +
                         snps_PC1 + snps_PC2 + snps_PC3 + snps_PC4 + snps_PC5)
fit <- lmFit(Mvalues_filt, design)
fit <- eBayes(fit)

DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)

pdf(file = "output_CrossPhase/Figures/QQplot_Baselineseverity_10CtrlPCs_5AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output_CrossPhase/Figures/QQplot_Baselineseverity_10CtrlPCs_5AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[7,2] <- lambda(DMPs$P.Value)

### 15 Ctrl PCs, 5 Anc PCs ###
design <- model.matrix(Mvalues_filt[1, ] ~ severity + Sex + Age + smoking + Collection_Year + cell_PC1 + cell_PC2 +
                         Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC5 + Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10 +
                         Ctrl_PC11 + Ctrl_PC12 + Ctrl_PC13 + Ctrl_PC14 + Ctrl_PC15 +
                         snps_PC1 + snps_PC2 + snps_PC3 + snps_PC4 + snps_PC5)
fit <- lmFit(Mvalues_filt, design)
fit <- eBayes(fit)

DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)

pdf(file = "output_CrossPhase/Figures/QQplot_Baselineseverity_15CtrlPCs_5AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output_CrossPhase/Figures/QQplot_Baselineseverity_15CtrlPCs_5AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[8,2] <- lambda(DMPs$P.Value)






### 2 Ctrl PCs, 10 Anc PCs ###
design <- model.matrix(Mvalues_filt[1, ] ~ severity + Sex + Age + smoking + Collection_Year + cell_PC1 + cell_PC2 + 
                         Ctrl_PC1 + Ctrl_PC2 + 
                         snps_PC1 + snps_PC2 + snps_PC3 + snps_PC4 + snps_PC5 + snps_PC6 + snps_PC7 + snps_PC8 + snps_PC9 + snps_PC10)
fit <- lmFit(Mvalues_filt, design)
fit <- eBayes(fit)

DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)

pdf(file = "output_CrossPhase/Figures/QQplot_Baselineseverity_2CtrlPCs_10AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output_CrossPhase/Figures/QQplot_Baselineseverity_2CtrlPCs_10AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[9,2] <- lambda(DMPs$P.Value)


### 5 Ctrl PCs, 10 Anc PCs ###
design <- model.matrix(Mvalues_filt[1, ] ~ severity + Sex + Age + smoking + Collection_Year + cell_PC1 + cell_PC2 + 
                         Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC5 +
                         snps_PC1 + snps_PC2 + snps_PC3 + snps_PC4 + snps_PC5 + snps_PC6 + snps_PC7 + snps_PC8 + snps_PC9 + snps_PC10)
fit <- lmFit(Mvalues_filt, design)
fit <- eBayes(fit)

DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)

pdf(file = "output_CrossPhase/Figures/QQplot_Baselineseverity_5CtrlPCs_10AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output_CrossPhase/Figures/QQplot_Baselineseverity_5CtrlPCs_10AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[10,2] <- lambda(DMPs$P.Value)


### 10 Ctrl PCs, 10 Anc PCs ###
design <- model.matrix(Mvalues_filt[1, ] ~ severity + Sex + Age + smoking + Collection_Year + cell_PC1 + cell_PC2 + 
                         Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC5 + Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10 +
                         snps_PC1 + snps_PC2 + snps_PC3 + snps_PC4 + snps_PC5 + snps_PC6 + snps_PC7 + snps_PC8 + snps_PC9 + snps_PC10)
fit <- lmFit(Mvalues_filt, design)
fit <- eBayes(fit)

DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)

pdf(file = "output_CrossPhase/Figures/QQplot_Baselineseverity_10CtrlPCs_10AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output_CrossPhase/Figures/QQplot_Baselineseverity_10CtrlPCs_10AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[11,2] <- lambda(DMPs$P.Value)

### 15 Ctrl PCs, 10 Anc PCs ###
design <- model.matrix(Mvalues_filt[1, ] ~ severity + Sex + Age + smoking + Collection_Year + cell_PC1 + cell_PC2 + 
                         Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC5 + Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10 +
                         Ctrl_PC11 + Ctrl_PC12 + Ctrl_PC13 + Ctrl_PC14 + Ctrl_PC15 +
                         snps_PC1 + snps_PC2 + snps_PC3 + snps_PC4 + snps_PC5 + snps_PC6 + snps_PC7 + snps_PC8 + snps_PC9 + snps_PC10)
fit <- lmFit(Mvalues_filt, design)
fit <- eBayes(fit)

DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)

pdf(file = "output_CrossPhase/Figures/QQplot_Baselineseverity_15CtrlPCs_10AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output_CrossPhase/Figures/QQplot_Baselineseverity_15CtrlPCs_10AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[12,2] <- lambda(DMPs$P.Value)




### 2 Ctrl PCs, 1 Anc PCs ###
design <- model.matrix(Mvalues_filt[1, ] ~ severity + Sex + Age + smoking + Collection_Year + cell_PC1 + cell_PC2 + 
                         Ctrl_PC1 + Ctrl_PC2 +
                         snps_PC1)
fit <- lmFit(Mvalues_filt, design)
fit <- eBayes(fit)

DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)

pdf(file = "output_CrossPhase/Figures/QQplot_Baselineseverity_2CtrlPCs_1AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output_CrossPhase/Figures/QQplot_Baselineseverity_2CtrlPCs_1AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[13,2] <- lambda(DMPs$P.Value)


### 5 Ctrl PCs, 1 Anc PCs ###
design <- model.matrix(Mvalues_filt[1, ] ~ severity + Sex + Age + smoking + Collection_Year + cell_PC1 + cell_PC2 +
                         Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC5 +
                         snps_PC1)
fit <- lmFit(Mvalues_filt, design)
fit <- eBayes(fit)

DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)

pdf(file = "output_CrossPhase/Figures/QQplot_Baselineseverity_5CtrlPCs_1AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output_CrossPhase/Figures/QQplot_Baselineseverity_5CtrlPCs_1AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[14,2] <- lambda(DMPs$P.Value)


### 10 Ctrl PCs, 1 Anc PCs ###
design <- model.matrix(Mvalues_filt[1, ] ~ severity + Sex + Age + smoking + Collection_Year + cell_PC1 + cell_PC2 +
                         Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC5 + Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10 +
                         snps_PC1)
fit <- lmFit(Mvalues_filt, design)
fit <- eBayes(fit)

DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)

pdf(file = "output_CrossPhase/Figures/QQplot_Baselineseverity_10CtrlPCs_1AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output_CrossPhase/Figures/QQplot_Baselineseverity_10CtrlPCs_1AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[15,2] <- lambda(DMPs$P.Value)


### 15 Ctrl PCs, 1 Anc PCs ###
design <- model.matrix(Mvalues_filt[1, ] ~ severity + Sex + Age + smoking + Collection_Year + cell_PC1 + cell_PC2 + 
                         Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC5 + Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10 +
                         Ctrl_PC11 + Ctrl_PC12 + Ctrl_PC13 + Ctrl_PC14 + Ctrl_PC15 +
                         snps_PC1)
fit <- lmFit(Mvalues_filt, design)
fit <- eBayes(fit)

DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)

pdf(file = "output_CrossPhase/Figures/QQplot_Baselineseverity_15CtrlPCs_1AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output_CrossPhase/Figures/QQplot_Baselineseverity_15CtrlPCs_1AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[16,2] <- lambda(DMPs$P.Value)



# EXPORT LAMBDA TABLE
write_xlsx(lambda_table, path = "output_CrossPhase/Tables/Lambda_Baselineseverity_limma.xlsx")



################################################################################
# annotate FINAL MODEL
################################################################################

### 5 Ctrl PCs, 1 Anc PC ###
design <- model.matrix(Mvalues_filt[1, ] ~ severity + Sex + Age + smoking + Collection_Year + 
                         cell_PC1 + cell_PC2 + Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC5 + snps_PC1)
fit <- lmFit(Mvalues_filt, design)
fit <- eBayes(fit)

DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
DMPs$cgID <- rownames(DMPs)
lambda(DMPs$P.Value) #0.9486756

# standard error
DMPs <- DMPs[order(rownames(DMPs)), ]
SE <- as.data.frame(sqrt(fit$s2.post) * fit$stdev.unscaled)
SE <- SE[order(rownames(SE)), ]
DMPs$SE <- SE$severity
DMPs <- DMPs[, c("cgID", "logFC", "CI.L", "CI.R", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "SE")]

# CORRECT FOR MULTIPLE TESTING
DMPs$adj_p_BH <- p.adjust(DMPs$P.Value, method = "BH")
DMPs$adj_p_Bonferroni <- p.adjust(DMPs$P.Value, method = "bonferroni")

# EXPORT
save(DMPs, file = "output_CrossPhase/RData/DMPs_Baselineseverity_5CtrlPCs_1AncPC_limma.RData")

# annotate RESULTS AND SAVE
DMPs_anno <- merge(DMPs, anno_upd, by.x = "cgID", by.y = "Name", all.x = TRUE)
save(DMPs_anno, file = "output_CrossPhase/RData/DMPs_Baselineseverity_5CtrlPCs_1AncPC_limma_anno.RData")
rm(DMPs)

# EXPORT SIGNIFICANT RESULTS
signif_DMPs_Bonferroni <- DMPs_anno[which(DMPs_anno$adj_p_Bonferroni <= 0.05),]
nrow(signif_DMPs_Bonferroni) #0
signif_DMPs_BH <- DMPs_anno[which(DMPs_anno$adj_p_BH <= 0.05),]
nrow(signif_DMPs_BH) #0
write_xlsx(signif_DMPs_Bonferroni, path = "output_CrossPhase/Tables/DMPs_Baselineseverity_5CtrlPCs_1AncPC_limma_signif_Bonferroni.xlsx")
save(signif_DMPs_Bonferroni, file = "output_CrossPhase/RData/DMPs_Baselineseverity_5CtrlPCs_1AncPC_limma_signif_Bonferroni.RData")
write_xlsx(signif_DMPs_BH, path = "output_CrossPhase/Tables/DMPs_Baselineseverity_5CtrlPCs_1AncPC_limma_signif_BH.xlsx")
save(signif_DMPs_BH, file = "output_CrossPhase/RData/DMPs_Baselineseverity_5CtrlPCs_1AncPC_limma_signif_BH.RData")




################################################################################
# DMRs
################################################################################

####################### PREPARE ################################################


stop("export summary statistics on Mac computer to run combp with ENmix package and import again")
# run:
# combp_input <- data.frame(chr = DMPs_anno$CHR, start = DMPs_anno$MAPINFO, end = DMPs_anno$MAPINFO + 1, p = DMPs_annp$P.Value, probe = DMPs_anno$cgID)
# combp_result <- combp(combp_input, dist.cutoff = 750, seed = 0.001, region_plot = F, mht_plot = F)


combp_result <- as.data.frame(read_xlsx("output_CrossPhase/Tables/combp_Severity_CrossPhase.xlsx"))

combp_sign <- adapt_combp(combp_result, DMPs_anno, anno_upd, nr_CpGs = 3, nr_signif_CpGs = 3, signif_proportion = 0.5)
#"significant DMRs before filtering step: 16"
#"significant DMRs after filtering: 5"

write_xlsx(combp_sign, "output_CrossPhase/Tables/DMRs_combp_CaseCtrl_5CtrlPCs_1AncPC_filtered.xlsx")


################################################################################
# MANHATTAN PLOT
################################################################################

# MANHATTAN PLOT
Manhattan_plot <- Manhattan(DMPs_anno, sex = "FM_auto", correction= "BH", combp_sign) #ignore warning messages

ggsave("output_CrossPhase/Figures/Manhattan_Baselineseverity_5CtrlPCs_1AncPC.svg", Manhattan_plot, units = "cm", width = 13, height = 9, dpi = 300)
ggsave("output_CrossPhase/Figures/Manhattan_Baselineseverity_5CtrlPCs_1AncPC.pdf", Manhattan_plot, units = "cm", width = 13, height = 9, dpi = 300)
ggsave("output_CrossPhase/Figures/Manhattan_Baselineseverity_5CtrlPCs_1AncPC.png", Manhattan_plot, units = "cm", width = 13, height = 9, dpi = 400)

