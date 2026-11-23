
setwd("/Users/kirahoeffler/OneDrive - University of Bergen/PhD/19_Meta-Analyses/Meta_Follow_Up/Severity_Blood/")

load("MetaResult_Blood_ProbesRemoved_Anno2.RData")

library(ENmix)

combp_input <- data.frame(chr = res.meta.anno$chr, start = res.meta.anno$pos, end = res.meta.anno$pos + 1, 
                          p = res.meta.anno$All_P_Random, probe = res.meta.anno$cgID)

combp_result <- combp(combp_input, dist.cutoff = 750, seed = 0.05, region_plot = F, mht_plot = F)
