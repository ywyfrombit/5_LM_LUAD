rm(list = ls())
library(rms)
library(pec)

set6 =  c("LYPD3","CCT6A","HMGA1","PSMB7","KRT8")
set1 = c("SH3BP5_AS1","LINC00654","MIR200CHG","TDRKH_AS1","TMPO_AS1","LINC00623","LINC01116",
         "ARNTL2_AS1","PCBP1_AS1","SNHG10","LINC00852","LINC00996","PLAC4","LINC01138")
set2 = c("PFKP","SLC2A1","BCAN","CDKN3","ANLN")
set3 = c("CPAMD8","HMMR","FSCN1","PKP2","KRT6A") 
set4 = c("AHSA1","SERBP1","RHOF","CCL20","CD3D")
set5 =  c("CACNA2D2", "CYP2B7P", "KRT6A")

library(GSVA)
library(survival)
library(DT)
library(survminer)
library(stringr)
library(dplyr)
library(ROCR)
library(caret)
library(ggplot2)
library(car)

setwd("G:/LUAD_scRNASeq/")


LUAD_TPM_clin = read.csv("G:/MR_mGWAS/LUAD_RNASeq/LUAD_TPM_Clin.csv",header = T,
                         row.names = 1,check.names = F)
LUAD_count = read.csv("G:/MR_mGWAS/LUAD_RNASeq/TCGA_LUAD_count_symbol.csv",row.names = 1,check.names = F)
LUAD_TPM = t(LUAD_TPM_clin[,12:ncol(LUAD_TPM_clin)])
colnames(LUAD_TPM) = LUAD_TPM_clin$X_PATIENT
load("cc_Wu_Zhou_2021_LMn.RData")
LMn = intersect(LMn$gene,rownames(LUAD_TPM))
LMn = list(LM = LMn)
LUAD_TPM_clin_fliter = LUAD_TPM_clin[!(LUAD_TPM_clin$OS ==0 & LUAD_TPM_clin$OS.time<30),]
LUAD_TPM_clin_fliter = LUAD_TPM_clin_fliter[!is.na(LUAD_TPM_clin_fliter$OS.time),]

library(ROCR)
library(caret)
library(ggplot2)
library(survminer)
library(car)

res_TCGA = data.frame(matrix(0,nrow=6,ncol=2))
colnames(res_TCGA) = c("Cindex","se")

for(i in 1:6){
  gene = get(paste0("set",as.character(i)))
  gene = intersect(gene,colnames(LUAD_TPM_clin_fliter))
  formula <- as.formula(paste("Surv(OS.time,OS) ~", paste(gene, collapse = " + ")))
  cox_fit <- coxph(formula, data = LUAD_TPM_clin_fliter)
  res_TCGA[i,1]=summary(cox_fit)$concordance[1]
  res_TCGA[i,2]=summary(cox_fit)$concordance[2]
}

setwd("E:/Bioinfor_reproduction/Rloop/microarray/LUAD")
load("./GSE13213/exprSet.Rdata")
load("./GSE13213/pd.Rdata")

exprSet = exprSet[,rownames(pd)[pd$Histology=="TP"]]
pd = pd[pd$Histology=="TP",]
clin_exp = data.frame(cbind(pd,t(exprSet)))

res_13213 = data.frame(matrix(0,nrow=6,ncol=2))
colnames(res_13213) = c("Cindex","se")
for(i in 1:6){
  gene = get(paste0("set",as.character(i)))
  gene = intersect(gene,colnames(clin_exp))
  formula <- as.formula(paste("Surv(OS.Time,OS) ~", paste(gene, collapse = " + ")))
  cox_fit <- coxph(formula, data = clin_exp)
  res_13213[i,1]=summary(cox_fit)$concordance[1]
  res_13213[i,2]=summary(cox_fit)$concordance[2]
}


setwd("E:/Bioinfor_reproduction/Rloop/microarray/LUAD")
load("./GSE26939/exprSet.Rdata")
load("./GSE26939/pd.Rdata")

exprSet = exprSet[,rownames(pd)[pd$Histology=="TP"]]
pd = pd[pd$Histology=="TP",]
clin_exp = data.frame(cbind(pd,t(exprSet)))
res_26939 = data.frame(matrix(0,nrow=6,ncol=2))
colnames(res_26939) = c("Cindex","se")
for(i in 2:6){
  gene = get(paste0("set",as.character(i)))
  gene = intersect(gene,colnames(clin_exp))
  formula <- as.formula(paste("Surv(OS.Time,OS) ~", paste(gene, collapse = " + ")))
  cox_fit <- coxph(formula, data = clin_exp)
  res_26939[i,1]=summary(cox_fit)$concordance[1]
  res_26939[i,2]=summary(cox_fit)$concordance[2]
}


setwd("E:/Bioinfor_reproduction/Rloop/microarray/LUAD")
load("./GSE30219/exprSet.Rdata")
load("./GSE30219/pd.Rdata")

exprSet = exprSet[,rownames(pd)[pd$Histology=="TP"]]
pd = pd[pd$Histology=="TP",]
clin_exp = data.frame(cbind(pd,t(exprSet)))
res_30219 = data.frame(matrix(0,nrow=6,ncol=2))
colnames(res_30219) = c("Cindex","se")
for(i in 1:6){
  gene = get(paste0("set",as.character(i)))
  gene = intersect(gene,colnames(clin_exp))
  formula <- as.formula(paste("Surv(OS.Time,OS) ~", paste(gene, collapse = " + ")))
  cox_fit <- coxph(formula, data = clin_exp)
  res_30219[i,1]=summary(cox_fit)$concordance[1]
  res_30219[i,2]=summary(cox_fit)$concordance[2]
}

setwd("E:/Bioinfor_reproduction/Rloop/microarray/LUAD")
load("./GSE31210/exprSet.Rdata")
load("./GSE31210/pd.Rdata")

exprSet = exprSet[,rownames(pd)[pd$Histology=="TP"]]
pd = pd[pd$Histology=="TP",]
clin_exp = data.frame(cbind(pd,t(exprSet)))
res_31210 = data.frame(matrix(0,nrow=6,ncol=2))
colnames(res_31210) = c("Cindex","se")
for(i in 1:6){
  gene = get(paste0("set",as.character(i)))
  gene = intersect(gene,colnames(clin_exp))
  formula <- as.formula(paste("Surv(OS.Time,OS) ~", paste(gene, collapse = " + ")))
  cox_fit <- coxph(formula, data = clin_exp)
  res_31210[i,1]=summary(cox_fit)$concordance[1]
  res_31210[i,2]=summary(cox_fit)$concordance[2]
}

setwd("E:/Bioinfor_reproduction/Rloop/microarray/LUAD")
load("./GSE41271/exprSet.Rdata")
load("./GSE41271/pd.Rdata")

exprSet = exprSet[,rownames(pd)[pd$Histology=="TP"]]
pd = pd[pd$Histology=="TP",]
clin_exp = data.frame(cbind(pd,t(exprSet)))
res_41271 = data.frame(matrix(0,nrow=6,ncol=2))
colnames(res_41271) = c("Cindex","se")
for(i in 1:6){
  gene = get(paste0("set",as.character(i)))
  gene = intersect(gene,colnames(clin_exp))
  formula <- as.formula(paste("Surv(OS.Time,OS) ~", paste(gene, collapse = " + ")))
  cox_fit <- coxph(formula, data = clin_exp)
  res_41271[i,1]=summary(cox_fit)$concordance[1]
  res_41271[i,2]=summary(cox_fit)$concordance[2]
}

res = cbind(res_TCGA$Cindex,res_13213$Cindex,res_26939$Cindex,res_30219$Cindex,
            res_31210$Cindex,res_41271$Cindex)

colnames(res) =c("TCGA-LUAD","GSE13213","GSE26939","GSE30219","GSE31210","GSE41271")
rownames(res) = c("Mai.S et al 2022", "Sun.J et al 2024","Zhao.F et al 2022",
                  "Zhang.p et al 2023","Guo.Z et al 2023","Ours")

rank_matrix <- apply(res, 2, function(x) rank(-x))
mean_ranks <- round(rowMeans(rank_matrix),2)
res_all <- cbind(res, Average_rank = mean_ranks)

write.csv(res_all,"G:/LUAD_scRNASeq/Fig_and_chart/LM_sig_LUAD_compare.csv",row.names = T)
