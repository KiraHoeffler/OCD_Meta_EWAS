
###
#DMRs identified using ENmix::combp on local machine (does not run on windows)
###


################################################################################
### SET UP ####
################################################################################

setwd("S:/Project/WP-epigenetics/13_MetaAnalyses/BaselineSeverity//")

#Libraries
library(writexl)
#library(dplyr)
#library(ggplot2)

load("output/RData/MetaResult_Blood_ProbesRemoved_Anno.RData")

#define candidate genes
cand_genes <- c("BDNF", "FKBP5", "NR3C1","OXTR", "SLC6A4") #MAOA -> did not look at sex chr

################################################################################
### FINALIZE DMRS ####
################################################################################

res.gene_comb <- as.data.frame(matrix(nrow = 0, ncol = 10))
colnames(res.gene_comb) <- c("cgID", "All_Effect_Random", "All_Effect_SE_Random", "All_P_Random", "P_Random_FDR", "chr", "pos", "UCSC_RefGene_Group", "padjust_gene", "Cand_Gene")

for (gene in cand_genes){
  res.gene <- res.meta.anno[grepl(gene, res.meta.anno$UCSC_RefGene_Name), ]
  res.gene$padjust_gene <- p.adjust(res.gene$All_P_Random)
  res.gene <- res.gene[, c("cgID", "All_Effect_Random", "All_Effect_SE_Random", "All_P_Random", "P_Random_FDR", "chr", "pos", "UCSC_RefGene_Group", "padjust_gene")]
  res.gene$Cand_Gene <- gene
  res.gene_comb <- rbind(res.gene_comb, res.gene)
}

write_xlsx(res.gene_comb, path = "output/Tables/Candidate_genes_blood.xlsx")

