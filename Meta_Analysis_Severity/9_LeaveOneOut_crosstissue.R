#################################################################################################
# Purpose of script:
################################################################################################

setwd("S:/Project/WP-epigenetics/13_MetaAnalyses/BaselineSeverity/")

library(doParallel)
library(foreach)
library(meta)
library(readxl)
library(purrr)
library(dplyr)
library(ggplot2)
library(forcats)
library(tidyr)
library(writexl)



load("output/RData/theme.RData")
load("output/RData/theme_transparent.Rdata")

#####################################################################
# Loading Data
#####################################################################

load("output/RData/MetaResult_CrossTissue_saliva_blood_ProbesRemoved_Anno_compr.RData")
meta_woVancouver <- res.meta.anno
load("input/UPDATED_ANNOTATION_v1v2.RData")

load("output/RData/MetaResult_Blood_ProbesRemoved_Anno_compr.RData")
blood <- res.meta.anno
blood <- blood[which(blood$cgID %in%  meta_woVancouver$cgID), ]
rm(res.meta.anno)
blood <- blood[, c("cgID", "All_Effect_Random", "All_Effect_SE_Random")]
colnames(blood) <- c("cgID", "logFC", "SE")
rownames(blood) <- blood$cgID

load("input/Bergen_DMPs_Baselineseverity_5CtrlPCs_1AncPC_limma.RData")
saliva <- DMPs
saliva <- saliva[which(saliva$cgID %in%  meta_woVancouver$cgID), ]
rm(DMPs)

#####################################################################
# Preparing 
#####################################################################
print("Setting up parallell")
cl = makeCluster(8)
registerDoParallel(cl)

#####################################################################
# SELECT PROBES PER DMR
#####################################################################



DMRs <- as.data.frame(read_xlsx("output/Tables/DMRs_combp_Meta_Severity_CrossTissue_saliva_blood.xlsx"))

CpG_list <- list()

for (i in c(1:nrow(DMRs))){
  chr <- DMRs$chr[i]
  start <- DMRs$start[i]
  end <- DMRs$end[i]
  DMR_cpgs <- unique(annov1v2$Name[which(annov1v2$chr == chr & annov1v2$pos >= start & annov1v2$pos <= end)])
  CpG_list[[i]] <- DMR_cpgs
}


names(CpG_list) <- DMRs$GeneString



sample_list <- list(full_meta = c("blood", "saliva"),
                    wo_blood = c("saliva"), 
                    wo_saliva = c("blood"))

res <- list(blood = blood, saliva = saliva)

#####################################################################
# Meta analysis
#####################################################################


print("Running the parallell loop for all probes to be analysed")
res.meta.DMRs = foreach(i=1:length(CpG_list), .combine=rbind, .packages = c("dplyr", "purrr")) %dopar% {
  DMR_name <- names(CpG_list)[i]
  probeIDs <- CpG_list[[i]]
  
  result_comb <- matrix(nrow = 0, ncol = 14)
  
  
  for (k in 1:length(sample_list)){
    
    result <- matrix(nrow = length(probeIDs), ncol = 14)
    
    result[, 1] <- names(sample_list)[k]
    result[, 2] <- DMR_name 
    result[, 3] <- probeIDs
    
    res_red <- res[sample_list[[k]]]
    
    
    for (j in 1:length(probeIDs)){
      probeID <- probeIDs[j]
      
      meanDiff <- sapply(res_red, function(x) x[["logFC"]][match(probeID,rownames(x))]) #logFC
      seDiff <- sapply(res_red, function(x) x[["SE"]][match(probeID,rownames(x))]) #SE
      
      meanDiff <- unlist(compact(meanDiff))
      seDiff <- unlist(compact(seDiff))
      
      out<- meta::metagen(TE = meanDiff, seTE = seDiff, method.tau = "REML")
      result[j, 4:14] <- c(sum(!is.na(meanDiff)), out$TE.fixed,out$seTE.fixed,out$pval.fixed, out$TE.random, out$seTE.random,out$pval.random, out$tau, out$I2, out$Q,1-pchisq(out$Q, out$df.Q))
      
    }
    
    result_comb <- rbind(result_comb, result)
  }
  
  return(result_comb)
}

## setting rownames and column names of result matrix
res.meta.DMRs <- as.data.frame(res.meta.DMRs)
colnames(res.meta.DMRs)<-c("Sample", "Gene", "cgID" , "N_cohorts", "All_Effect_Fixed", "All_Effect_SE_Fixed", "All_P_Fixed", "All_Effect_Random", "All_Effect_SE_Random","All_P_Random", "All_tau", "All_I2", "All_Q", "All_Het P")
res.meta.DMRs <- res.meta.DMRs[which(res.meta.DMRs$N_cohorts != 0), ]

res.meta.DMRs <- res.meta.DMRs %>%
  mutate(across(-c(Gene, Sample, cgID), ~ if(is.character(.x)) suppressWarnings(as.numeric(.x)) else .x
  ))

# Saving result
save(res.meta.DMRs,file="output/RData/MetaResult_CrossTissue_DMRs_LeaveOneOut_woVancouver")




plot_df <- res.meta.DMRs %>%
  group_by(Gene, Sample) %>%
  summarise(
    n_cpg = n_distinct (cgID),
    estimate = median(All_Effect_Random),
    low = quantile(All_Effect_Random, 0.25), 
    high = quantile(All_Effect_Random, 0.75),
    frac_positive = mean(All_Effect_Random > 0),
    .groups = "drop")

gene_sizes <- plot_df %>%
  group_by(Gene) %>%
  summarise(n_cpg_total = max(n_cpg), .groups = "drop")

plot_df <- plot_df %>%
  left_join(gene_sizes, by = "Gene") %>%
  mutate(
    facet_label = paste0(Gene, " (", n_cpg_total, " CpGs)")
  )

effect_plot <- ggplot(plot_df, aes(x = estimate, y = fct_rev(Sample))) + 
  geom_vline(xintercept = 0, linewidth = 0.5) + 
  geom_errorbarh(aes(xmin = low, xmax = high, height = 0.2)) + 
  geom_point(size = 2) + 
  facet_wrap(~facet_label, ncol = 3, scales = "free_y") + 
  labs(
    x = "median CpG meta-analysis effect size (± IQR)",
    y = NULL, 
    size = "CpGs", 
    title = "Leave-one-out meta analysis"
  ) + th + th_transparent + 
  theme(strip.text = element_text(face = "bold", size = 10))  + 
  theme(strip.text = element_text(hjust = 0.5, face = "bold"))

effect_plot

ggsave(plot = effect_plot, filename =  "output/Figures/Effects_Meta_Severity_CrossTissue_woVancouver.png", width = 22, height = 6, units = "cm", dpi = 500)
ggsave(plot = effect_plot, filename =  "output/Figures/Effects_Meta_Severity_CrossTissue_woVancouver.pdf", width = 22, height = 6, units = "cm")



sign_df <- res.meta.DMRs %>%
  mutate(nominal_sign = All_P_Random < 0.05) %>%
  group_by(Gene, Sample) %>%
  summarise(
    n_cpg = n_distinct(cgID),
    n_nominal = sum(nominal_sign),
    frac_nominal = 100 * (n_nominal/n_cpg), 
    .groups = "drop"
  )

gene_sizes <- sign_df %>%
  group_by(Gene) %>%
  summarise(n_cpg_total = max(n_cpg), .groups = "drop")

sign_df <- sign_df %>%
  left_join(gene_sizes, by="Gene") %>%
  mutate(
    facet_label = paste0(Gene, " (", n_cpg_total, " CpGs)")
  )


LOO_CpG_plot <- ggplot(sign_df, aes(x=frac_nominal, y=fct_rev(Sample))) + 
  geom_col(width = 0.7) + 
  facet_wrap(~facet_label, ncol = 2, scales = "free_y") +
  labs(x = "% nominally significant CpGs (p < 0.05)", y=NULL, title = "Nominally significant CpGs in LOO") + 
  th + th_transparent +
  theme(strip.text = element_text(face = "bold", size = 10)) + 
  theme(strip.text = element_text(hjust = 0.5, face = "bold"))
LOO_CpG_plot

ggsave(plot = LOO_CpG_plot, filename =  "output/Figures/SignificantCpGs_Meta_Severity_CrossTissue_woVancouver.png", width = 17, height = 8, units = "cm", dpi = 500)
ggsave(plot = LOO_CpG_plot, filename =  "output/Figures/SignificantCpGs_Meta_Severity_CrossTissue_woVancouver.pdf", width = 17, height = 8, units = "cm")


library(tidyr)

influence_df <- plot_df %>%
  select(Gene, Sample, estimate) %>%
  mutate(estimate = abs(estimate)) %>%
  pivot_wider(names_from = Sample, values_from = estimate) %>%
  mutate(across(starts_with("wo"), 
                ~ .x - full_meta, 
                .names = "delta_{.col}")) %>%
  pivot_longer(
    starts_with("delta_"),
    names_to = "LOO",
    values_to = "delta"
  ) %>%
  mutate(
    LOO = sub("^delta_", "", LOO)
  )

influence_df <- influence_df %>%
  mutate(
    Gene = factor(Gene, levels = rev(sort(unique(Gene))))
  )

Effect_change_plot <- ggplot(influence_df, aes(x=LOO, y= Gene, fill = delta)) + 
  geom_tile() +
  ggtitle("Change in DMR effect size in LOO") + 
  scale_fill_gradient2(low = "#25abd0" , mid = "white", high =  "#d04125") + th + th_transparent +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank()) + 
  theme(strip.text = element_text(hjust = 0.5, face = "bold"))
Effect_change_plot

ggsave(plot = Effect_change_plot, filename =  "output/Figures/Effect_Change_Meta_Severity_CrossTissue_woVancouver.png", width = 12, height = 10, units = "cm", dpi = 500)
ggsave(plot = Effect_change_plot, filename =  "output/Figures/Effect_Change_Meta_Severity_CrossTissue_woVancouver.pdf", width = 12, height = 10, units = "cm")






res_long <- imap_dfr(res, ~ .x %>% mutate(cohort = .y))

res_DMR_cohorts <- foreach(i = 1:length(CpG_list), .combine = rbind, .packages = c("dplyr", "purrr")) %dopar%{
  DMR_name <- names(CpG_list)[i]
  probeIDs <- CpG_list[[i]]
  n_CpGs <- length(probeIDs)
  
  res_long %>%
    filter(cgID %in% probeIDs) %>%
    group_by(cohort) %>%
    summarise(
      DMR = paste0(DMR_name, " (", n_CpGs, " CpGs)"),
      median_effect = median(logFC), q1 = quantile(logFC, 0.25), q3 = quantile(logFC, 0.75), .groups = "drop"
    )
}


plot_df2 <- plot_df[which(plot_df$Sample == "full_meta"), ]
plot_df2 <- plot_df2[, c( "Sample", "facet_label", "estimate", "low", "high")]
colnames(plot_df2) <- colnames(res_DMR_cohorts)
plot_df2$DMR <- c("NUDT12 (8 CpGs)", "PEMT (7 CpGs)", "TTC23 (8 CpGs)")
res_DMR_cohorts2 <- rbind(res_DMR_cohorts, plot_df2)
res_DMR_cohorts2$cohort[which(res_DMR_cohorts2$cohort == "full_meta")] <- "meta"
res_DMR_cohorts2$cohort <- factor(res_DMR_cohorts2$cohort, levels = c( "meta", "saliva", "blood"))


Effect_cohort_plot <- ggplot(res_DMR_cohorts2, aes(x= median_effect, y= cohort)) + 
  geom_vline(xintercept = 0, linewidth = 0.4) + 
  geom_point(size = 2) + 
  geom_errorbarh(aes(xmin = q1, xmax = q3), height = 0.2) + 
  facet_wrap(~DMR, scales = "free_x", ncol = 2) + 
  labs(
    title = "Cohort-specific DMR effect sizes", 
    x = "median CpG effect size per DMR (± IQR)",
    y = NULL
  ) + th + th_transparent + 
  theme(strip.text = element_text(hjust = 0.5, face = "bold"))

Effect_cohort_plot

ggsave(plot = Effect_cohort_plot, filename =  "output/Figures/Effect_Cohort_Meta_Severity_CrossTissue_woVancouver.png", width = 13, height = 8.5, units = "cm", dpi = 500)
ggsave(plot = Effect_cohort_plot, filename =  "output/Figures/Effect_Cohort_Meta_Severity_CrossTissue_woVancouver.pdf", width = 13, height = 8.5, units = "cm")

