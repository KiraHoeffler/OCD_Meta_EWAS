
################################################################################
# LIBRARIES #
################################################################################

library(missMethyl)
library(writexl)
library(ggplot2)

setwd("/Users/kirahoeffler/OneDrive - University of Bergen/PhD/19_Meta-Analyses/Meta_Follow_Up/Severity_Blood/")

################################################################################
# LOAD & ADAPT ANNOTATION
################################################################################

load("../UPDATED_ANNOTATION_v1v2.RData")
anno_upd <- unique(annov1v2[, c("Name", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "chr", "pos")])
colnames(anno_upd) <- c("Name", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "chr", "pos")
anno_upd <- anno_upd[which(anno_upd$chr != "chr0"), ]
rownames(anno_upd) <- anno_upd$Name

################################################################################
# LOAD AND ADAPT SUMMARY STATS FOR GO #
################################################################################

load("MetaResult_Blood_ProbesRemoved_Anno_compr.RData")

all_cpgs <- res.meta.anno$cgID
ordered_DMPs_anno <- res.meta.anno[order(res.meta.anno$All_P_Random, decreasing = FALSE), ]
sign1percent <- round(nrow(res.meta.anno)*0.01, 0)

signDMPs <- ordered_DMPs_anno$cgID[1:sign1percent]
#sign_0005 <- ordered_DMPs_anno$cgID[which(ordered_DMPs_anno$P.Value < 0.005)]

################################################################################
# RUN GO AND ADAPT RESULTS #
################################################################################

gometh <- function (sig.cpg, all.cpg = NULL, collection = c("GO", "KEGG"), 
          array.type = c("450K", "EPIC", "EPIC_V2"), plot.bias = FALSE, 
          prior.prob = TRUE, anno = NULL, equiv.cpg = TRUE, fract.counts = TRUE, 
          genomic.features = c("ALL", "TSS200", "TSS1500", "Body", 
                               "1stExon", "3'UTR", "5'UTR", "ExonBnd"), sig.genes = FALSE) 
{
  array.type <- match.arg(toupper(array.type), c("450K", "EPIC", 
                                                 "EPIC_V2"))
  collection <- match.arg(toupper(collection), c("GO", "KEGG"))
  genomic.features <- match.arg(genomic.features, c("ALL", 
                                                    "TSS200", "TSS1500", "Body", "1stExon", "3'UTR", "5'UTR", 
                                                    "ExonBnd"), several.ok = TRUE)
  if (length(genomic.features) > 1 & any(grepl("ALL", genomic.features))) {
    message("All input CpGs are used for testing.")
    genomic.features <- "ALL"
  }
  if (array.type == "450K" & any(grepl("ExonBnd", genomic.features))) {
    stop("'ExonBnd' is not an annotated feature on 450K arrays,\n\n           please remove it from your genomic.feature parameter\n\n           specification.")
  }
  if (array.type == "EPIC_V2" & any(grepl("ExonBnd", genomic.features))) {
    stop("'ExonBnd' is not an annotated feature on EPIC_V2 arrays,\n\n           please remove it from your genomic.feature parameter\n\n           specification.")
  }
  if (collection == "GO") {
    
      if(!requireNamespace("org.Hs.eg.db", quietly = TRUE))
        stop("org.Hs.eg.db package required but not installed.")
      egGO2ALLEGS <- utils::getFromNamespace("org.Hs.egGO2ALLEGS", "org.Hs.eg.db")
      GeneID.PathID <- AnnotationDbi::toTable(egGO2ALLEGS)[,c("gene_id", "go_id", "Ontology")]
      d <- !duplicated(GeneID.PathID[, c("gene_id", "go_id")])
      GeneID.PathID <- GeneID.PathID[d, ]
      GOID.TERM <- suppressMessages(AnnotationDbi::select(GO.db::GO.db, 
                                                          keys=unique(GeneID.PathID$go_id), 
                                                          columns=c("GOID","ONTOLOGY","TERM"), 
                                                          keytype="GOID"))
      go <- tapply(GeneID.PathID$gene_id, GeneID.PathID$go_id, list)
      
      go <- list(idList=go, idTable=GOID.TERM)
      
      go$idTable <- go$idTable[which(go$idTable$ONTOLOGY == "BP"), ]
      
      go$idList <- go$idList[names(go$idList) %in% go$idTable$GOID]
      
    result <- gsameth(sig.cpg = sig.cpg, all.cpg = all.cpg, 
                      collection = go$idList, array.type = array.type, 
                      plot.bias = plot.bias, prior.prob = prior.prob, anno = anno, 
                      equiv.cpg = equiv.cpg, fract.counts = fract.counts, 
                      genomic.features = genomic.features, sig.genes = sig.genes)
    result <- merge(go$idTable, result, by.x = "GOID", by.y = "row.names")
    rownames(result) <- result$GOID
  }
  else if (collection == "KEGG") {
    kegg <- .getKEGG()
    result <- gsameth(sig.cpg = sig.cpg, all.cpg = all.cpg, 
                      collection = kegg$idList, array.type = array.type, 
                      plot.bias = plot.bias, prior.prob = prior.prob, anno = anno, 
                      equiv.cpg = equiv.cpg, fract.counts = fract.counts, 
                      genomic.features = genomic.features, sig.genes = sig.genes)
    result <- merge(kegg$idTable, result, by.x = "PathwayID", 
                    by.y = "row.names")
    rownames(result) <- result$PathwayID
  }
  result[, -1]
}



GO_result <- gometh(sig.cpg = signDMPs, all.cpg = all_cpgs, collection = "GO", array.type = "EPIC", genomic.features = "ALL", anno = anno_upd)
GO_result$logP <- -log10(GO_result$P.DE)
GO_result <- GO_result[order(GO_result$P.DE), ]


################################################################################
# VISUALIZE & EXPORT RESULTS #
################################################################################

# top 5 results
results <- GO_result[which(GO_result$N <3000), ]
results$significant <- ifelse(results$FDR <= 0.05, "Y", "N")

write_xlsx(results, path = "/Users/kirahoeffler/OneDrive - University of Bergen/PhD/19_Meta-Analyses/Meta_Follow_Up/Severity_Blood/GO.xlsx")

results_red2 <- results[1:10, ]

GO_plot <- ggplot(results_red2, aes(x = logP, y = reorder(TERM, logP), fill = significant)) +
  geom_bar(stat = "identity") +
  labs(x = "-logP", y = NULL, title = NULL) +
  ggtitle("GO Enrichment for the Top 1% Differentially Methylated CpGs") +
  scale_fill_manual(values = c("Y" =  "#192f2d","N" = "#B0D0CE")) +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(0, NA), expand = c(0,0)) + 
  theme(plot.title = element_text(size = 14, color = "black", face = "bold"),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 0, size = 12, vjust = 0.5, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 12, vjust = 0.5, colour = "black"),
        axis.title.x = element_text(vjust = -0.5, size = 12, colour = "black"),
        axis.title.y = element_text(vjust = -0.5, size = 12, colour = "black", margin = margin(r = 10)),
        axis.title = element_text(face = "bold"),
        panel.border = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size =12),
        legend.key = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent"))
GO_plot

# export
ggsave("/Users/kirahoeffler/OneDrive - University of Bergen/PhD/19_Meta-Analyses/Meta_Follow_Up//Severity_Blood/GO_Meta_Severity_Blood.pdf", plot = GO_plot, device = "pdf", width = 22, height = 11, units = "cm", dpi = 400)
