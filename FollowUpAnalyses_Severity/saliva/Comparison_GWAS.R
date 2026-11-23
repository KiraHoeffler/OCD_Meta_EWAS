
library(readxl)


set.seed(123)

# Libraries
library(regioneR)
library(readxl)
library(writexl)

setwd(dir = "/Users/kirahoeffler/OneDrive - University of Bergen/PhD/19_Meta-Analyses/Meta_Follow_Up/Severity_Saliva/")

Strom <- read.table("/Users/kirahoeffler/OneDrive - University of Bergen/PhD/14_Downloads/GWAS_OCS/26352157/ocs2024obsessive-compulsive_symptoms_daner_STR_NTR_SfS_TwinsUK_strometal", header = TRUE)
load("../UPDATED_ANNOTATION_v1v2.RData")

annov1v2$chr[which(annov1v2$chr == "chr14_GL000009v2_random")] <- "chr14"
annov1v2$chr[which(annov1v2$chr == "chr14_KI270846v1_alt")] <- "chr14"
annov1v2$chr[which(annov1v2$chr == "chr17_KI270857v1_alt")] <- "chr17"
annov1v2$chr[which(annov1v2$chr == "chr17_KI270909v1_alt")] <- "chr17"
annov1v2$chr[which(annov1v2$chr == "chr19_KI270938v1_alt")] <- "chr19"
annov1v2$chr[which(annov1v2$chr == "chr22_KI270879v1_alt")] <- "chr22"
annov1v2$chr[which(annov1v2$chr == "chr4_GL000008v2_random")] <- "chr4"
annov1v2$chr[which(annov1v2$chr == "chr8_KI270821v1_alt")] <- "chr8"


### LOWEST P VALUE ###
DMRs_Meta <- as.data.frame(read_xlsx("DMRs_Baselineseverity_5CtrlPCs_1AncPC_extended.xlsx"))
DMRs_Meta$GWAS_lowestP <- NA

for (i in 1:nrow(DMRs_Meta)){
  chr <- as.numeric(gsub("chr", "", DMRs_Meta$chr[i]))
  start_pos <- DMRs_Meta$start[i]
  end_pos <- DMRs_Meta$end[i]
  
  Strom_red <- Strom[which(Strom$CHR == chr & (Strom$BP > start_pos - 10000) & (Strom$BP < end_pos + 10000) ), ]
  DMRs_Meta$GWAS_lowestP[i] <- min(Strom_red$P)
}




# PERMUTATION ###


DMRs_Meta$GWAS_lowestP <- NA
DMRs_Meta$GWAS_highest_logOR <- NA

for (i in 1:nrow(DMRs_Meta)){
  chr <- as.numeric(gsub("chr", "", DMRs_Meta$chr[i]))
  start_pos <- DMRs_Meta$start[i]
  end_pos <- DMRs_Meta$end[i]
  
  Strom_red <- Strom[which(Strom$CHR == chr & (Strom$BP > start_pos - 1000) & (Strom$BP < end_pos + 1000) ), ]
  Strom_red$abslogOR <- abs(log(Strom_red$OR))
  
  if (nrow(Strom_red) > 0){
    DMRs_Meta$GWAS_lowestP[i] <- min(Strom_red$P)
    DMRs_Meta$GWAS_highest_logOR[i] <- max(Strom_red$abslogOR)
  } else {
    DMRs_Meta$GWAS_lowestP[i] <- NA
  }
}



genome <- getGenomeAndMask(data.frame(CHR = paste0("chr", c(1:22, "X", "Y")), 
                                      chr_length = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067,
                                                     159138663, 146364022, 141213431, 135534747, 135006516, 133851895,
                                                     115169878, 107349540, 102531392, 90354753, 81195210, 78077248,
                                                     59128983, 63025520, 48129895, 51304566, 0, 0)))


DMRs_Meta$p_GWAS_permutation <- NA
DMRs_Meta$z_GWAS_permutation <- NA

for (j in 1:nrow(DMRs_Meta)){
  
  if(!is.na(DMRs_Meta$GWAS_lowestP[j])){
    print(j)
    # Convert to GRanges
    DMRs_gr <- toGRanges(DMRs_Meta[j, c("chr", "start", "end")])
    
    random_list <- replicate(
      1500,
      randomizeRegions(
        A = DMRs_gr,
        genome = genome$genome, 
        mask = DMRs_gr
      ),
      simplify = FALSE
    )
    
    rand_pval <- c()
    for (i in 1:length(random_list)){
      regions <- as.data.frame(random_list[[i]], row.names = NULL)
      
      chr <- as.numeric(gsub("chr", "", regions$seqnames[1]))
      start_pos <- regions$start[1]
      end_pos <- regions$end[1]
      
      Strom_red <- Strom[which(Strom$CHR == chr & (Strom$BP > start_pos - 1000) & (Strom$BP < end_pos + 1000) ), ]
      
      if (nrow(Strom_red) > 0){
        rand_pval <- c(rand_pval, min(Strom_red$P))
      } 
    }
    
    
    real_pval <- -log10(DMRs_Meta$GWAS_lowestP[j])
    rand_pval <- -log10(rand_pval[!is.na(rand_pval)])
    
    p_val_perm <- (sum(rand_pval >= real_pval)) / (length(rand_pval))
    z_score <- (real_pval - mean(rand_pval))/sd(rand_pval)
    
    DMRs_Meta$p_GWAS_permutation[j] <- p_val_perm
    DMRs_Meta$z_GWAS_permutation[j] <- z_score
  }
}

write_xlsx(DMRs_Meta, path = "Comparison_GWAS.xlsx")




