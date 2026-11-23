
library(mashr)
library(dplyr)
library(readxl)
library(writexl)
library(ggplot2)
library(tidyr)

setwd("S:/Project/WP-epigenetics/13_MetaAnalyses/BaselineSeverity//")

load("output/RData/theme.RData")

####################################################################################
# Read DMP Summary Stats
####################################################################################

load("output/RData/MetaResult_Blood_ProbesRemoved_Anno.RData")
blood <- res.meta.anno[, c("cgID", "All_Effect_Random", "All_Effect_SE_Random", "All_P_Random")]
colnames(blood) <- c("cgID", "Effect_blood", "SE_blood", "P.Value_blood")

load("input/Bergen_DMPs_Baselineseverity_5CtrlPCs_1AncPC_limma.RData")
saliva <- DMPs[, c("cgID", "logFC", "SE", "P.Value")]
colnames(saliva) <- c("cgID", "Effect_saliva", "SE_saliva", "P.Value_saliva")

rm(res.meta.anno, DMPs)

# Merge
merged_df <- merge(blood, saliva, by= "cgID")
merged_df <- na.omit(merged_df)

####################################################################################
# Read DMRs
####################################################################################

blood_DMRs <- as.data.frame(read_xlsx("output/Tables/DMRs_DMPs_blood.xlsx"))
saliva_DMRs <- as.data.frame(read_xlsx("../../11_BaselineSeverity/output_CrossPhase/Tables/DMRs_Baselineseverity_5CtrlPCs_1AncPC_extended.xlsx"))

load("input/UPDATED_ANNOTATION_v1v2.RData")

blood_DMRs$median_p_saliva <- NA
blood_DMRs$median_effect <- NA
blood_DMRs$nom_sign_percentage <- NA
blood_DMRs$effect_direction_percentage <- NA

for (i in 1:nrow(blood_DMRs)){
  chr <- blood_DMRs$chr[i]
  start <- blood_DMRs$start[i]
  end <- blood_DMRs$end[i]
  
  DMR_cpgs <- unique(annov1v2$Name[which(annov1v2$chr == chr & annov1v2$pos >= start & annov1v2$pos <= end)])
  
  merged_df_red <- merged_df[which(merged_df$cgID %in% DMR_cpgs), ]
  
  blood_DMRs$median_p_saliva[i] <- median(merged_df_red$P.Value_saliva)
  blood_DMRs$median_effect[i] <- median(merged_df_red$Effect_saliva)
  blood_DMRs$nom_sign_percentage[i] <- (nrow(merged_df_red[which(merged_df_red$P.Value_saliva < 0.05), ]))/nrow(merged_df_red)
  blood_DMRs$effect_direction_percentage[i] <- (length(merged_df_red %>% filter(sign(Effect_blood) == sign(Effect_saliva)) %>% pull(cgID)))/(nrow(merged_df_red))
}


saliva_DMRs$median_p_blood <- NA
saliva_DMRs$median_effect <- NA
saliva_DMRs$nom_sign_percentage <- NA
saliva_DMRs$effect_direction_percentage <- NA

for (i in 1:nrow(saliva_DMRs)){
  chr <- saliva_DMRs$chr[i]
  start <- saliva_DMRs$start[i]
  end <- saliva_DMRs$end[i]
  
  DMR_cpgs <- unique(annov1v2$Name[which(annov1v2$chr == chr & annov1v2$pos >= start & annov1v2$pos <= end)])
  
  merged_df_red <- merged_df[which(merged_df$cgID %in% DMR_cpgs), ]
  
  saliva_DMRs$median_p_blood[i] <- median(merged_df_red$P.Value_blood)
  saliva_DMRs$median_effect[i] <- median(merged_df_red$Effect_blood)
  saliva_DMRs$nom_sign_percentage[i] <- (nrow(merged_df_red[which(merged_df_red$P.Value_blood < 0.05), ]))/nrow(merged_df_red)
  saliva_DMRs$effect_direction_percentage[i] <- (length(merged_df_red %>% filter(sign(Effect_saliva) == sign(Effect_blood)) %>% pull(cgID)))/(nrow(merged_df_red))
}


write_xlsx(saliva_DMRs, path = "output/Tables/Concordance_saliva_DMRs.xlsx")
write_xlsx(blood_DMRs, path = "output/Tables/Concordance_blood_DMRs.xlsx")


cor(merged_df$Effect_blood, merged_df$Effect_saliva, method = "spearman") #0.01585675
mean(sign(merged_df$Effect_blood) == sign(merged_df$Effect_saliva)) #0.5077819


plot_df_saliva <- saliva_DMRs %>%
  mutate(locus = factor(GeneString, levels = saliva_DMRs$GeneString[order(saliva_DMRs$effect_direction_percentage, decreasing = TRUE)]),
         nominally_significant = nom_sign_percentage * 100,
         same_effect_direction = effect_direction_percentage *100) %>%
  select(locus, nominally_significant, same_effect_direction) %>%
  pivot_longer(
    cols = c(nominally_significant, same_effect_direction), 
    names_to = "category", values_to = "percentage"
  )

saliva_plot <- ggplot(plot_df_saliva, aes(x=factor(locus), y = percentage, fill = category)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = c(nominally_significant = "#E69F00", same_effect_direction = "#0072B2"), 
                    labels = c(nominally_significant = "nominally significant CpGs in blood", 
                               same_effect_direction = "CpGs with same direction of effects in blood")) +
  ggtitle("Saliva DMRs - concordance with blood") + 
  labs(x = "locus", y = "percentage", fill = NULL) + 
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100,20), labels = function(x) paste0(x, "%")) + 
  th + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position ="top", legend.title = element_blank()) + 
  guides(fill = guide_legend(ncol = 1))
saliva_plot  


ggsave(plot = saliva_plot, filename = "output/Figures/saliva_concordance.pdf", width = 13, height = 13, units = "cm")






plot_df_blood <- blood_DMRs %>%
  mutate(locus = factor(GeneString, levels = blood_DMRs$GeneString[order(blood_DMRs$effect_direction_percentage, decreasing = TRUE)]),
         nominally_significant = nom_sign_percentage * 100,
         same_effect_direction = effect_direction_percentage *100) %>%
  select(locus, nominally_significant, same_effect_direction) %>%
  pivot_longer(
    cols = c(nominally_significant, same_effect_direction), 
    names_to = "category", values_to = "percentage"
  )

blood_plot <- ggplot(plot_df_blood, aes(x=factor(locus), y = percentage, fill = category)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = c(nominally_significant = "#E69F00", same_effect_direction = "#0072B2"), 
                    labels = c(nominally_significant = "nominally significant CpGs in saliva", 
                               same_effect_direction = "CpGs with same direction of effects in saliva")) +
  ggtitle("blood DMRs/DMP - concordance with saliva") + 
  labs(x = "locus", y = "percentage", fill = NULL) + 
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100,20), labels = function(x) paste0(x, "%")) + 
  th + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position ="top", legend.title = element_blank()) + 
  guides(fill = guide_legend(ncol = 1))
blood_plot  


ggsave(plot = blood_plot, filename = "output/Figures/blood_concordance.pdf", width = 13, height = 13, units = "cm")


merged_df_sign_saliva <- merged_df[which(merged_df$P.Value_saliva < 0.05), ]
nrow(merged_df_sign_saliva[which(merged_df_sign_saliva$P.Value_blood < 0.05), ])/nrow(merged_df_sign_saliva) #0.04593153

merged_df_sign_blood <- merged_df[which(merged_df$P.Value_blood < 0.05), ]
nrow(merged_df_sign_blood[which(merged_df_sign_blood$P.Value_saliva < 0.05), ])/nrow(merged_df_sign_blood) #0.04705177
