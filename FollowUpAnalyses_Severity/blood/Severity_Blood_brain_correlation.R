
set.seed(123)

# Libraries
library(regioneR)
library(readxl)
library(writexl)
library(dplyr)

setwd("/Users/kirahoeffler/OneDrive - University of Bergen/PhD/19_Meta-Analyses/Meta_Follow_Up/Severity_Blood//")

#####################################################################################

load("/Users/kirahoeffler/OneDrive - University of Bergen/PhD/19_Meta-Analyses/Meta_Follow_Up/UPDATED_ANNOTATION_v1v2.RData")
our_DMRs <- as.data.frame(read_xlsx("DMRs_DMPs.xlsx"))

con <- gzfile("MetaResult_Blood_ProbesRemoved_Anno_compr.RData", "rb")
load(con)
close(con)


annov1v2 <- annov1v2[which(annov1v2$Name %in% res.meta.anno$cgID), ]
annov1v2$chr[which(annov1v2$chr == "chr14_GL000009v2_random")] <- "chr14"
annov1v2$chr[which(annov1v2$chr == "chr14_KI270846v1_alt")] <- "chr14"
annov1v2$chr[which(annov1v2$chr == "chr17_KI270857v1_alt")] <- "chr17"
annov1v2$chr[which(annov1v2$chr == "chr17_KI270909v1_alt")] <- "chr17"
annov1v2$chr[which(annov1v2$chr == "chr19_KI270938v1_alt")] <- "chr19"
annov1v2$chr[which(annov1v2$chr == "chr22_KI270879v1_alt")] <- "chr22"
annov1v2$chr[which(annov1v2$chr == "chr4_GL000008v2_random")] <- "chr4"
annov1v2$chr[which(annov1v2$chr == "chr8_KI270821v1_alt")] <- "chr8"

# sort CpGs
cpg_annot <- as.data.frame(annov1v2 %>%
                             arrange(chr, pos))
cpg_annot$pos <- as.numeric(cpg_annot$pos)

# distance between CpGs
cpg_annot$distance <- ave(
  cpg_annot$pos,
  cpg_annot$chr,
  FUN = function(x) c(NA, diff(x))
)

# start a new region when first CpG in chromosome or >750bp to previous CpG, give regions numbers
cpg_annot$region_id <- ave(
  cpg_annot$distance,
  cpg_annot$chr,
  FUN = function(x) cumsum(is.na(x) | x > 750)
)

#summarise nr of cpgs per region
region_summary <- cpg_annot %>%
  group_by(chr, region_id) %>%
  summarise(
    start = min(pos), 
    end = max(pos),
    n_cpgs = n(),
    .groups = "drop"
  )

dmrs <- as.data.frame(region_summary %>% filter(n_cpgs >= 3))
nrow(dmrs)
dmrs <- dmrs %>% arrange(chr, start)
save(dmrs, file = "/Users/kirahoeffler/OneDrive - University of Bergen/PhD/19_Meta-Analyses/Meta_Follow_Up/Possible_DMRs_blood.RData")


chrom_sizes <- data.frame(chr = paste0("chr", c(1:22, "X", "Y")), 
                         chr_length = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067,
                                        159138663, 146364022, 141213431, 135534747, 135006516, 133851895,
                                        115169878, 107349540, 102531392, 90354753, 81195210, 78077248,
                                        59128983, 63025520, 48129895, 51304566, 155270560, 59373566))

# function to identify regions not overlapping with DMRs
find_non_dmrs <- function(chr_dmrs, chr_length){
 #if no DMRs on chr -> select whole chr
 if (nrow(chr_dmrs) == 0){
   return(data.frame(start = 1, end = chr_length))
 }
 non_regions <- list()
 
 #region before first DMR
 first_start <- chr_dmrs$start[1]
 if(first_start > 1){
   region_before <- data.frame(start = 1, end = first_start - 1)
   non_regions[[length(non_regions) + 1]] <- region_before
 }
 
 #gaps between DMRs
 for (i in 1:(nrow(chr_dmrs) - 1)){
   gap_start <- chr_dmrs$end[i] + 1
   gap_end <- chr_dmrs$start[i + 1] - 1
   
   if(gap_start <= gap_end){
     region_between <- data.frame(start = gap_start, end = gap_end)
     non_regions[[length(non_regions) + 1]] <- region_between
   }
 }
 
 #region after last DMR
 last_end <- chr_dmrs$end[nrow(chr_dmrs)]
 if (last_end < chr_length) {
   region_after <- data.frame(start = last_end + 1, end = chr_length)
   non_regions[[length(non_regions) + 1]] <- region_after
 }
 
 #combine
 non_dmrs_chr <- do.call(rbind, non_regions)
 return(non_dmrs_chr)
}


#apply to the whole annotation
dmrs <- na.omit(dmrs)

non_dmrs <- as.data.frame(dmrs %>%
 group_by(chr) %>%
 group_map(~ {
   chr <- .y$chr
   chr_length <- chrom_sizes$chr_length[chrom_sizes$chr == chr]
   chr_dmrs <- .x
   non_regions <- find_non_dmrs(chr_dmrs, chr_length)
   non_regions$chr <- chr
   non_regions
 }) %>%
 bind_rows() %>%
 select(chr, start, end))

save(non_dmrs, file = "../NonDMRs_Blood.RData")




#####################################################################################
# Import
corr_table <- read.table("/Users/kirahoeffler/Library/CloudStorage/OneDrive-UniversityofBergen/PhD/14_Downloads/ImageCpG/Illimina_EPIC_covar", header = TRUE, sep = " ", stringsAsFactors = FALSE)

genome <- getGenomeAndMask(data.frame(chr = paste0("chr", c(1:22)), 
                         chr_length = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067,
                         159138663, 146364022, 141213431, 135534747, 135006516, 133851895,
                         115169878, 107349540, 102531392, 90354753, 81195210, 78077248,
                         59128983, 63025520, 48129895, 51304566)))

our_DMRs$brain_rho <- NA
our_DMRs$brain_p <- NA
our_DMRs$brain_sign_any <- NA
our_DMRs$brain_median_corr <- NA
for (i in 1:nrow(our_DMRs)){
  cpgs <- unlist(strsplit(our_DMRs$probe[i], split = ";"))
  extr_cgs <- corr_table[which(corr_table$cgid %in% cpgs), ]
  sign_any <- any(extr_cgs$p_rho_brain_saliva <= 0.05)
  
  our_DMRs$brain_rho[i] <- paste0(round(extr_cgs$rho_brain_saliva, 3), collapse = ";")
  our_DMRs$brain_p[i] <- paste0(round(extr_cgs$p_rho_brain_saliva, 3), collapse = ";")
  our_DMRs$brain_sign_any[i] <- sign_any
  our_DMRs$brain_median_corr[i] <- abs(median(extr_cgs$rho_brain_saliva, na.rm = TRUE))
}

our_DMRs$brain_rho_blood <- NA
our_DMRs$brain_p_blood <- NA
our_DMRs$brain_sign_any_blood <- NA
our_DMRs$brain_median_corr_blood <- NA
for (i in 1:nrow(our_DMRs)){
  cpgs <- unlist(strsplit(our_DMRs$probe[i], split = ";"))
  extr_cgs <- corr_table[which(corr_table$cgid %in% cpgs), ]
  sign_any <- any(extr_cgs$p_rho_brain_blood <= 0.05)
  
  our_DMRs$brain_rho_blood[i] <- paste0(round(extr_cgs$rho_brain_blood, 3), collapse = ";")
  our_DMRs$brain_p_blood[i] <- paste0(round(extr_cgs$p_rho_brain_blood, 3), collapse = ";")
  our_DMRs$brain_sign_any_blood[i] <- sign_any
  our_DMRs$brain_median_corr_blood[i] <- abs(median(extr_cgs$rho_brain_blood, na.rm = TRUE))
}
                         
non_dmrs <- non_dmrs[which(non_dmrs$chr != "chrX" & non_dmrs$chr != "chrY"),]
                         
our_DMRs$p_saliva_brain_permutation <- NA
our_DMRs$z_saliva_brain_permutation <- NA

our_DMRs$p_blood_brain_permutation <- NA
our_DMRs$z_blood_brain_permutation <- NA

for (j in 1:nrow(our_DMRs)){
  print(j)
  # Convert to GRanges
  DMRs_gr <- toGRanges(our_DMRs[j, c("chr", "start", "end")])
  noDMR_gr <- toGRanges(non_dmrs[, c("chr", "start", "end")])
  
  random_list <- replicate(
    1500,
    randomizeRegions(
      A = DMRs_gr,
      genome = genome$genome, 
      mask = noDMR_gr
    ),
    simplify = FALSE
  )
  
  rand_median <- c()
  for (i in 1:length(random_list)){
    regions <- as.data.frame(random_list[[i]], row.names = NULL)
    
    cpgs <- annov1v2$Name[which(annov1v2$chr == as.character(regions$seqnames[1]) & 
                                  annov1v2$pos >= regions$start[1] & 
                                  annov1v2$pos <= regions$end[1])]
    extr_cgs <- corr_table[which(corr_table$cgid %in% cpgs), ]
    brain_median_corr <- abs(median(extr_cgs$rho_brain_saliva, na.rm = TRUE))
    rand_median <- c(rand_median, brain_median_corr)
  }
  
  
  rand_median_blood <- c()
  for (i in 1:length(random_list)){
    regions <- as.data.frame(random_list[[i]], row.names = NULL)
    
    cpgs <- annov1v2$Name[which(annov1v2$chr == as.character(regions$seqnames[1]) & 
                                  annov1v2$pos >= regions$start[1] & 
                                  annov1v2$pos <= regions$end[1])]
    extr_cgs <- corr_table[which(corr_table$cgid %in% cpgs), ]
    brain_median_corr <- abs(median(extr_cgs$rho_brain_blood, na.rm = TRUE))
    rand_median_blood <- c(rand_median_blood, brain_median_corr)
  }
  
  
  real_median <- our_DMRs$brain_median_corr[j]
  rand_median <- rand_median[!is.na(rand_median)]
  
  real_median_blood <- our_DMRs$brain_median_corr_blood[j]
  rand_median_blood <- rand_median_blood[!is.na(rand_median_blood)]
  
  p_val_perm <- (sum(rand_median >= real_median)) / (length(rand_median))
  z_score <- (real_median - mean(rand_median))/sd(rand_median)
  
  p_val_perm_blood <- (sum(rand_median_blood >= real_median_blood)) / (length(rand_median_blood))
  z_score_blood <- (real_median_blood - mean(rand_median_blood))/sd(rand_median_blood)
  
  our_DMRs$p_saliva_brain_permutation[j] <- p_val_perm
  our_DMRs$z_saliva_brain_permutation[j] <- z_score
  
  our_DMRs$p_blood_brain_permutation[j] <- p_val_perm_blood
  our_DMRs$z_blood_brain_permutation[j] <- z_score_blood
  
}

write_xlsx(our_DMRs, path = "/Users/kirahoeffler/OneDrive - University of Bergen/PhD/19_Meta-Analyses/Meta_Follow_Up/Severity_Blood/Severity_blood_correlation_brain.xlsx")



########### DMPs ###############

# Import
load("MetaResult_Blood_ProbesRemoved_Anno_compr.RData")
DMPs <- res.meta.anno[which(res.meta.anno$P_Random_FDR < 0.05), ]

extr_DMPs <- corr_table[which(corr_table$cgid %in% DMPs$cgID), ]

DMPs$brain_rho <- NA
DMPs$brain_p <- NA
DMPs$brain_sign_any <- NA
for (i in 1:nrow(DMPs)){
  cpgs <- DMPs$cgID[i]
  extr_cgs <- corr_table[which(corr_table$cgid %in% cpgs), ]
  sign_any <- any(extr_cgs$p_rho_brain_blood <= 0.05)
  
  DMPs$brain_rho[i] <- paste0(round(extr_cgs$rho_brain_blood, 3), collapse = ";")
  DMPs$brain_p[i] <- paste0(round(extr_cgs$p_rho_brain_blood, 3), collapse = ";")
  DMPs$brain_sign_any[i] <- sign_any
}

DMPs$p_blood_brain_permutation <- NA
DMPs$z_blood_brain_permutation <- NA

for(i in 1:nrow(DMPs)){
  print(i)
  cgs_minus_ours <- setdiff(annov1v2$Name, DMPs$cgID[i])
  
  random_cpgs <- sample(cgs_minus_ours, 1500)
  
  real_corr <- abs(as.numeric(DMPs$brain_rho[i]))
  random_corr <- abs(corr_table$rho_brain_blood[which(corr_table$cgid %in% random_cpgs)])
  random_corr <- random_corr[!is.na(random_corr)]
  
  
  p_val_perm <- (sum(random_corr >= real_corr) + 1) / (length(random_corr) + 1)
  z_score <- (real_corr - mean(random_corr))/sd(random_corr)
  
  DMPs$p_blood_brain_permutation[i] <- p_val_perm
  DMPs$z_blood_brain_permutation[i] <- z_score
}

write_xlsx(DMPs, path = "/Users/kirahoeffler/OneDrive - University of Bergen/PhD/19_Meta-Analyses/Meta_Follow_Up/Severity_Blood/Severity_blood_correlation_brain_DMP.xlsx")
