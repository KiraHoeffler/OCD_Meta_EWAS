#################################################################################################
# Purpose of script: Meta-analysis. 
################################################################################################

setwd("S:/Project/WP-epigenetics/13_MetaAnalyses/BaselineSeverity/")

library(doParallel)
library(foreach)
library(meta)
library(tictoc)

#####################################################################
# Input
#####################################################################

#saliva
load("input/Bergen_DMPs_Baselineseverity_5CtrlPCs_1AncPC_limma.RData")
saliva <- DMPs[, c("cgID", "logFC", "SE")]
colnames(saliva) <- c("cgID", "Effect", "SE")
rownames(saliva) <- saliva$cgID


#blood
load("output/RData/MetaResult_Blood_ProbesRemoved_Anno.RData")
blood <- res.meta.anno[, c("cgID", "All_Effect_Random", "All_Effect_SE_Random")]
colnames(blood) <- c("cgID", "Effect", "SE")
rownames(blood) <- blood$cgID


#combine
res <- list(blood, saliva)
names(res) = c("blood", "saliva")
save(res, file = "output/RData/MetaInput_CrossTissue_saliva_blood.RData")

#####################################################################
# Preparing 
#####################################################################
print("Setting up parallell")
cl = makeCluster(10)
registerDoParallel(cl)

#####################################################################
# Meta analysis
#####################################################################
print("Setting up analysis")
n_cohorts = length(res)
probes = unlist(sapply(res,rownames))
probes<-table(probes)
probes<-probes[which(probes >= 2)]

tic()
print("Running the parallell loop for all probes to be analysed")
res.meta = foreach(i=1:length(probes), .combine=rbind) %dopar% {
  probeID<-names(probes)[i]
  
  meanDiff <- sapply(res, function(x) x[["Effect"]][match(probeID,rownames(x))]) #logFC
  seDiff <- sapply(res, function(x) x[["SE"]][match(probeID,rownames(x))]) #SE
  
  tryCatch({
    out<- meta::metagen(TE = meanDiff, seTE = seDiff, method.tau = "REML")
    return(c(sum(!is.na(meanDiff)), out$TE.fixed,out$seTE.fixed,out$pval.fixed, out$TE.random, out$seTE.random,out$pval.random, out$tau, out$I2, out$Q,1-pchisq(out$Q, out$df.Q)))
  },
  error = function(e) print(paste("error",i,sep=": ")))
}
toc() #7846.57 sec elapsed

## setting rownames and column names of result matrix
res.meta <- as.data.frame(res.meta)
rownames(res.meta)<-names(probes)
colnames(res.meta)<-c("N_cohorts", "All_Effect_Fixed", "All_Effect_SE_Fixed", "All_P_Fixed", "All_Effect_Random", "All_Effect_SE_Random","All_P_Random", "All_tau", "All_I2", "All_Q", "All_Het P")
res.meta$P_Random_FDR <- p.adjust(res.meta$All_P_Random)

# Saving result
print("Saving result")
save(res.meta,file="output/RData/MetaResult_CrossTissue_saliva_blood.RData")

print("Done")
