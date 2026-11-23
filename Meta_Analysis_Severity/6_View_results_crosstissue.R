#################################################################################################
# Purpose of script: Removing probes that should have been removed, creating QQ plots and saving top lists
################################################################################################

setwd("S:/Project/WP-epigenetics/13_MetaAnalyses/BaselineSeverity/")

library(qqman)
library(Cairo)

# LAMBDA FUNCTION
calculate_lambda <- function(pvalues){
  chisq <- qchisq(1-pvalues, 1)
  lambda <- median(chisq)/qchisq(0.5, 1)
  return(lambda)
}

######################################################################################################
# Check cohorts for probes that should have been removed
######################################################################################################
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
anno.epic2 = getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
anno.epic = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno.450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

load("../../02_Import/OCD_Pipeline/Resources/450K_RemoveProbes_ChenProbeIDs.Robj")
remove450 = remove_probes
load("../../02_Import/OCD_Pipeline/Resources/EPICv1_RemoveProbes.Robj")
removeEpic = remove_probes
load("../../02_Import/OCD_Pipeline/Resources/EPICv2_RemoveProbes.RData")
removeEpic2 = remove_probes
removeEpic2 = unique(gsub("_.*$", "", removeEpic2))

######################################################################################################
# View Results
######################################################################################################

load("output/RData/MetaResult_CrossTissue_saliva_blood.RData")
res.meta$All_P_Fixed = as.numeric(res.meta$All_P_Fixed)
res.meta$All_P_Random = as.numeric(res.meta$All_P_Random)

if(any(is.na(res.meta$All_P_Random)))
  res.meta = res.meta[-which(is.na(res.meta$All_P_Random)),]
if(length(intersect(rownames(res.meta),remove450))>0)
  res.meta = res.meta[-which(rownames(res.meta) %in% intersect(rownames(res.meta),remove450)),]
if(length(intersect(rownames(res.meta),removeEpic))>0)
  res.meta = res.meta[-which(rownames(res.meta) %in% intersect(rownames(res.meta),removeEpic)),]
if(length(intersect(rownames(res.meta),removeEpic2))>0)
  res.meta = res.meta[-which(rownames(res.meta) %in% intersect(rownames(res.meta),removeEpic2)),]

res.meta$P_Random_FDR <- p.adjust(res.meta$All_P_Random)

save(res.meta,file = "output/RData/MetaResult_CrossTissue_saliva_blood_ProbesRemoved.RData")


# ANNOTATE
load("input/UPDATED_ANNOTATION_v1v2.RData")
res.meta$cgID <- rownames(res.meta)
res.meta.anno <- merge(res.meta, annov1v2, by.x = "cgID", by.y = "Name", all.x = TRUE)

res.meta.anno$All_Effect_Fixed <- as.numeric(res.meta.anno$All_Effect_Fixed)
res.meta.anno$All_Effect_SE_Fixed <- as.numeric(res.meta.anno$All_Effect_SE_Fixed)
res.meta.anno$All_Effect_Random <- as.numeric(res.meta.anno$All_Effect_Random)
res.meta.anno$All_Effect_SE_Random <- as.numeric(res.meta.anno$All_Effect_SE_Random)
res.meta.anno$All_tau <- as.numeric(res.meta.anno$All_tau)
res.meta.anno$All_I2 <- as.numeric(res.meta.anno$All_I2)
res.meta.anno$All_Q <- as.numeric(res.meta.anno$All_Q)
res.meta.anno$`All_Het P` <- as.numeric(res.meta.anno$`All_Het P`)
res.meta.anno$pos <- as.numeric(res.meta.anno$pos)

save(res.meta.anno,file = "output/RData/MetaResult_CrossTissue_saliva_blood_ProbesRemoved_Anno.RData")
save(res.meta.anno,file = "output/RData/MetaResult_CrossTissue_saliva_blood_ProbesRemoved_Anno_compr.RData", version = 3, compress = "gzip")



### QQ plots ##
# Fixed
pdf(file = "output/Figures/QQplot_CrossTissue_saliva_blood_fixed.pdf", width = 5, height = 5)
qq(res.meta$All_P_Fixed)
dev.off()

calculate_lambda(res.meta$All_P_Fixed) #0.9730508

# Random
pdf(file = "output/Figures/QQplot_CrossTissue_saliva_blood_random.pdf", width = 5, height = 5)
qq(res.meta$All_P_Random)
dev.off()

png(file = "output/Figures/QQplot_CrossTissue_saliva_blood_random2.png", width = 5, height = 5, units = "in", res = 300)
qq(res.meta$All_P_Random)
dev.off()

calculate_lambda(res.meta$All_P_Random) #0.8022608



# Saving top CpGs with pval < 1e-6
res.meta.top = res.meta[which(res.meta$All_P_Random<1e-6),]
write.csv(res.meta.top,file = "output/Tables//Meta_crosstissue_saliva_blood_randomTOP1e-6.csv",row.names = T,quote = F)


