
###
#DMRs identified using ENmix::combp on local machine (does not run on windows)
###


################################################################################
### SET UP ####
################################################################################

setwd("S:/Project/WP-epigenetics/13_MetaAnalyses/BaselineSeverity//")

#Libraries
library(writexl)
library(dplyr)
library(ggplot2)

# Load
load("output/RData/theme.RData")
load("output/RData/theme_transparent.Rdata")

combp_result <- as.data.frame(read.csv("output/Tables/resu_combp_Severity_Blood.csv"))
load("output/RData/MetaResult_Blood_ProbesRemoved_Anno.RData")
load("input/UPDATED_ANNOTATION_v1v2.RData")

################################################################################
### FINALIZE DMRS ####
################################################################################

print(paste0("significant DMRs before filtering step: ",nrow(combp_result))) #47
  
# only select significant DMRs and those with at least 3 CpGs
combp_sign <- combp_result[which(as.numeric(combp_result$sidak) <= 0.05 & as.numeric(combp_result$nprobe) >= 3), ]
  
# nominal significant DMPs
res.meta.anno_sign05 <- res.meta.anno[which(res.meta.anno$All_P_Random < 0.05), ]
  
#annotate DMP_overlap
combp_sign$GeneString <- NA
combp_sign$RegulRegion <- NA
combp_sign$Nr_CpGs <- NA
combp_sign$Nr_signif_CpGs_05 <- NA
  
  
for (i in c(1:nrow(combp_sign))){
  chr <- combp_sign$chr[i]
  start <- combp_sign$start[i]
  end <- combp_sign$end[i]
  
  Genes <- unique(annov1v2$UCSC_RefGene_Name[which(annov1v2$chr == chr & annov1v2$pos >= start &  annov1v2$pos <= end)])
  Genes <- unique(unlist(strsplit(Genes, ";")))
  GeneString <- paste(Genes, collapse = ";")
  combp_sign$GeneString[i] <- GeneString
  
  RegulRegion <- unique(annov1v2$UCSC_RefGene_Group[which(annov1v2$chr == chr & annov1v2$pos >= start &  annov1v2$pos <= end)])
  RegulRegion <- unique(unlist(strsplit(RegulRegion, ";")))
  RegulRegionString <- paste(RegulRegion, collapse = ";")
  combp_sign$RegulRegion[i] <- RegulRegionString
  
  DMR_cpgs <- unique(annov1v2$Name[which(annov1v2$chr == chr & annov1v2$pos >= start & annov1v2$pos <= end)])
  combp_sign$Nr_CpGs[i] <- length(DMR_cpgs)
  nr_sign05 <- length(intersect(DMR_cpgs, res.meta.anno_sign05$cgID))
  combp_sign$Nr_signif_CpGs_05[i] <- nr_sign05
}
  
combp_sign <- combp_sign[which(combp_sign$Nr_signif_CpGs_05 >= 3), ]
  
combp_sign$Perc_sign_CpGs <- combp_sign$Nr_signif_CpGs_05 / (combp_sign$Nr_CpGs)
combp_sign <- combp_sign[which(combp_sign$Perc_sign_CpGs >= 0.5), ]
  
print(paste0("significant DMRs after filtering: ", nrow(combp_sign))) #10
  
# oder: adj. p value
combp_sign <- combp_sign[order(combp_sign$sidak), ]
  
# add UNANNOTATED
length_unannotated <- length(combp_sign$GeneString[which(combp_sign$GeneString == "")])
  
count <- 1
for (i in c(1:nrow(combp_sign))){
  if (combp_sign$GeneString[i] == ""){
    New_Gene <- paste0("UNANNOTATED", count)
    combp_sign$GeneString[i] <- New_Gene
    count <- count + 1
  }
}

write_xlsx(combp_sign, "output/Tables/DMRs_combp_Meta_Severity_Blood.xlsx")


################################################################################
### MANHATTAN PLOT ####
################################################################################


# PREPARATION
all_DMR_cpgs <- c()

for (i in c(1:nrow(combp_sign))){
  chr <- combp_sign$chr[i]
  start <- combp_sign$start[i]
  end <- combp_sign$end[i]
  DMR_cpgs <- unique(annov1v2$Name[which(annov1v2$chr == chr & annov1v2$pos >= start & annov1v2$pos <= end)])
  all_DMR_cpgs <- c(all_DMR_cpgs, DMR_cpgs)
}
rm (DMR_cpgs)

all_DMR_cpgs_pvalues <- res.meta.anno$All_P_Random[which(res.meta.anno$cgID %in% all_DMR_cpgs)]
min_DMR_p <- -log10(min(all_DMR_cpgs_pvalues))
rm(all_DMR_cpgs_pvalues)


# MANHATTAN PLOT
signif_DMPs <- res.meta.anno[which(res.meta.anno$P_Random_FDR <= 0.05),]
min_sign <- -log10(max(signif_DMPs$All_P_Random))

DMPs_anno <- res.meta.anno

# ADD -LOG10 P VALUES, MAPINFO, AND COLOR TO THE RESULTS
DMPs_anno$pvalue_Manhattan <- -log10(DMPs_anno$All_P_Random)
DMPs_anno$color <- NA
DMPs_anno$chr <- substr(DMPs_anno$chr, 4, nchar(DMPs_anno$chr))
odd <- c("1","3","5","7","9","11","13","15","17","19","21","X")
DMPs_anno$color <- ifelse(DMPs_anno$chr %in% odd, "#d4e8e6", "#a8d1ce")
DMPs_anno$color[which(DMPs_anno$cgID %in% all_DMR_cpgs)] <- "#3c726d"
DMPs_anno$color[which(DMPs_anno$cgID %in% signif_DMPs$cgID)] <- "#192f2d"
DMPs_anno$chr2 <- paste0("chr", DMPs_anno$chr)


DMPs_anno$chr2 <- factor(DMPs_anno$chr2, levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                                                    "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
                                                     "chr19", "chr20", "chr21", "chr22"))


#continuous counting of positions
count_data <- DMPs_anno %>%
  group_by(chr2) %>%
  summarise(max_position = max(pos)) %>%
  mutate(position_add = lag(cumsum(as.numeric(max_position)), default = 0)) %>%
  select(chr2, position_add)

#add actual position with continous chromosome counting position
DMPs_anno2 <- DMPs_anno %>%
  inner_join(count_data, by = "chr2") %>%
  mutate(position_cum = pos + position_add)

#position for x axis chr label
axis_point <- DMPs_anno2 %>%
  group_by(chr2) %>%
  summarise(center = mean(position_cum))

axis_point <- na.omit(axis_point)

axis_point$chr <- c("1","2","3","4","5","6","7","8","","10","","12","","14","","16","","18","","","","")

# y limit
y_limit <- max(DMPs_anno2$pvalue_Manhattan) + 0.5

# DMR coordinates
DMR_coordinates <- DMPs_anno2$position_cum[which(DMPs_anno2$color == "#3c726d")]

# ACTUAL PLOT
Manhattan_plot <- ggplot() +
  geom_point(data = subset(DMPs_anno2, color == "#d4e8e6"), aes(x=position_cum, y=pvalue_Manhattan, color = color), shape = 19, size = 0.8) +
  geom_point(data = subset(DMPs_anno2, color == "#a8d1ce"), aes(x=position_cum, y=pvalue_Manhattan, color = color), shape = 19, size = 0.8) +
  geom_point(data = subset(DMPs_anno2, color == "#3c726d"), aes(x=position_cum, y=pvalue_Manhattan, color = color), shape = 19, size = 0.8) +
  geom_point(data = subset(DMPs_anno2, color == "#192f2d"), aes(x=position_cum, y=pvalue_Manhattan, color = color), shape = 19, size = 0.8) +
  scale_size_manual(values = c(1,2)) +
  th + th_transparent + 
  theme(axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.x = element_text(size = 10)) +
  labs(x = "genomic position (chromosome)", y = "-log10(p)") +
  scale_x_continuous(expand = c(0,0), label = axis_point$chr, breaks = axis_point$center) + 
  scale_y_continuous(expand = c(0,0), limits = c(0, y_limit)) +
  scale_color_identity() +
  geom_hline(yintercept = 5, color = "darkgrey", linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = min_sign, color = "darkgrey", linetype = "dashed", alpha = 0.5) +
  lapply(DMR_coordinates, function(x) geom_segment(aes(x = x, xend = x, y = 0, yend = min_DMR_p), linewidth = 0.25, color = "#3c726d"))

ggsave("output/Figures/Manhattan_Meta_Severity_Blood.pdf", Manhattan_plot, units = "cm", width = 14, height = 9, dpi = 300)
ggsave("output/Figures/Manhattan_Meta_Severity_Blood.png", Manhattan_plot, units = "cm", width = 14, height = 9, dpi = 400)
