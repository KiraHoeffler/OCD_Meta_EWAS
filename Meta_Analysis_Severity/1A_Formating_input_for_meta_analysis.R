#################################################################################################
# Purpose of script: Load summary statistics from individual studies and put in a list formatted for meta-analysis
################################################################################################

setwd("S:/Project/WP-epigenetics/13_MetaAnalyses/BaselineSeverity/")

library(qqman)

################################################################################################
# SAMPLES
################################################################################################

load("input/Barcelona_DMPs_BaselineSeverity_Sex_adjusted_1CtrlPCs_1AncPCs_limma.RData")
Barcelona = DMPs
Barcelona$cgID <- rownames(DMPs)
load("input/Bergen_DMPs_Baselineseverity_5CtrlPCs_1AncPC_limma.RData")
Bergen = DMPs
Bergen$cgID <- rownames(DMPs)
load("input/Berlin_DMPs_BaselineSeverity_Sex_adjusted_5CtrlPCs_2AncPCs_limma.RData")
Berlin = DMPs
Berlin$cgID <- rownames(DMPs)
load("input/Bonn_DMPs_BaselineSeverity_Sex_adjusted_15CtrlPCs_2AncPCs_limma.RData")
Bonn = DMPs
Bonn$cgID <- rownames(DMPs)
load("input/Freiburg_DMPs_BaselineSeverity_Sex_adjusted_5CtrlPCs_2AncPCs_limma.RData")
Freiburg = DMPs
Freiburg$cgID <- rownames(DMPs)
load("input/Vancouver_DMPs_BaselineSeverity_Sex_adjusted_2CtrlPCs_2AncPCs_limma.RData")
Vancouver = DMPs
Vancouver$cgID <- rownames(DMPs)


################################################################################################
# Lambda
################################################################################################

# LAMBDA
source("input/lambda.R")
lambda(Barcelona$P.Value) #0.8807302
lambda(Bergen$P.Value) #0.9486756
lambda(Berlin$P.Value) #1.029668
lambda(Bonn$P.Value) #1.053372
lambda(Freiburg$P.Value) #1.061966
lambda(Vancouver$P.Value) #1.038251

# QQ PLOTS
pdf(file = "output//Figures/QQplot_Barcelona.pdf", width = 5, height = 5)
qq(Barcelona$P.Value)
dev.off()

pdf(file = "output//Figures/QQplot_Bergen.pdf", width = 5, height = 5)
qq(Bergen$P.Value)
dev.off()

pdf(file = "output//Figures/QQplot_Berlin.pdf", width = 5, height = 5)
qq(Berlin$P.Value)
dev.off()

pdf(file = "output//Figures/QQplot_Bonn.pdf", width = 5, height = 5)
qq(Bonn$P.Value)
dev.off()

pdf(file = "output//Figures/QQplot_Freiburg.pdf", width = 5, height = 5)
qq(Freiburg$P.Value)
dev.off()

pdf(file = "output//Figures/QQplot_Vancouver.pdf", width = 5, height = 5)
qq(Vancouver$P.Value)
dev.off()


################################################################################################
# Formatting for all samples
################################################################################################

res = list(Barcelona, Bergen, Berlin, Bonn, Freiburg, Vancouver)
names(res) = c("Barcelona", "Bergen", "Berlin", "Bonn", "Freiburg", "Vancouver")
save(res, file="output/RData/MetaInput_AllSamples.RData")


################################################################################################
# Formatting for blood only
################################################################################################

res = list(Barcelona, Berlin, Bonn, Freiburg)
names(res) = c("Barcelona", "Berlin", "Bonn", "Freiburg")
save(res, file="output/RData/MetaInput_BloodSamples.RData")


################################################################################################
# Formatting for leaving out Freiburg from Blood
################################################################################################

res = list(Barcelona, Berlin, Bonn)
names(res) = c("Barcelona", "Berlin", "Bonn")
save(res, file="output/RData/MetaInput_BloodSamples_woFreiburg.RData")


