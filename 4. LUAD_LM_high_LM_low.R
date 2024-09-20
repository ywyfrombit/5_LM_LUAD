library(GSVA)
library(survival)
library(DT)
library(survminer)
library(stringr)
library(dplyr)
setwd("G:/LUAD_scRNASeq/")


LUAD_TPM_clin = read.csv("G:/MR_mGWAS/LUAD_RNASeq/LUAD_TPM_Clin.csv",header = T,
                         row.names = 1,check.names = F)
LUAD_count = read.csv("G:/MR_mGWAS/LUAD_RNASeq/TCGA_LUAD_count_symbol.csv",row.names = 1,check.names = F)
LUAD_TPM = t(LUAD_TPM_clin[,12:ncol(LUAD_TPM_clin)])
colnames(LUAD_TPM) = LUAD_TPM_clin$X_PATIENT
load("cc_Wu_Zhou_2021_LMn.RData")
LMn = intersect(LMn$gene,rownames(LUAD_TPM))
LMn = list(LM = LMn)

#LMn GSVA
LM_gsva <- gsva(expr = LUAD_TPM, gset.idx.list = LMn, kcdf = "Gaussian",parallel.sz = 8)
LM_gsva = as.data.frame(t(LM_gsva))
LM_gsva$LM_group = ifelse(LM_gsva$LM > median(LM_gsva$LM),"high_LM","low_LM")
all(rownames(LM_gsva) == LUAD_TPM_clin$X_PATIENT)
LUAD_TPM_clin$LM_group = LM_gsva$LM_group


fit <- survfit(Surv(OS.time, OS) ~ LM_group, 
               data = LUAD_TPM_clin)
x=ggsurvplot(fit,
             pval = TRUE, 
             pval.coord = c(1430, 0.9),
             main = "Survival curve",
             conf.int = F,
             legend.title = "Cluster", legend.labs = c("LM_high","LM_low"), 
             break.time.by = 365,
             xlim = c(0, 3650),
             risk.table = TRUE, # Add risk table
             risk.table.col = "strata", # Change risk table color by groups
             tables.height = 0.2, 
             tables.theme = theme_cleantable(), 
             linetype = "strata", # Change line type by groups
             surv.median.line = "hv" # Specify median survival
             )
  
x$plot = x$plot + ggtitle("TCGA-LUAD") + theme(plot.title = element_text(size = 15, hjust = 0.5)) 
x

data.survdiff <- survdiff(Surv(OS.time,OS)~ LM_group ,data = LUAD_TPM_clin)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
ci <- paste0(sprintf("%.3f",HR)," [",sprintf("%.3f",low95),", ",sprintf("%.3f",up95),"]")

x$plot <- x$plot + 
  annotate("text", x = 1400, y = 0.83, 
           label = paste0("\n HR (95% CI) = ",ci),  
           size = 5, color = "black", hjust = 0)+
  theme(text = element_text(size = 15))
x$plot


########################################################################################################

#LDHA, KRT6A, VDAC1, GAPDH, LYPD3

library(ROCR)
library(caret)
library(ggplot2)
library(survminer)
library(car)

fit_B <- coxph(Surv(OS.time,OS)~ CCT6A+KRT8+LYPD3+HMGA1+PSMB7 ,data = LUAD_TPM_clin)
survminer::ggforest(fit_B)
vif_values = vif(fit_B)
barplot(vif_values, main = "VIF Values", horiz = TRUE, col = "steelblue")

library(timeROC)
risk_B = predict(fit_B,newx=LUAD_TPM_clin[,c("CCT6A","KRT8","LYPD3","HMGA1","PSMB7")],type="lp")

LUAD_TPM_clin$risk_B = risk_B
LUAD_TPM_clin$group_B = ifelse(LUAD_TPM_clin$risk_B>median(LUAD_TPM_clin$risk_B),"high-risk","low-risk")
surv_fit_B = survfit(Surv(OS.time,OS)~ group_B ,data = LUAD_TPM_clin)
x=ggsurvplot(surv_fit_B,
             pval = TRUE,
             pval.coord = c(1500, 0.9),
             main = "Survival curve",
             conf.int = F,
             break.time.by = 365,
             xlim = c(0, 3650),
             risk.table = TRUE, # Add risk table
             risk.table.col = "strata", # Change risk table color by groups
             tables.height = 0.2,
             tables.theme = theme_cleantable(), 
             linetype = "strata", # Change line type by groups
             surv.median.line = "hv" # Specify median survival
)
x$plot = x$plot + ggtitle("LM risk gene in TCGA-LUAD") + theme(plot.title = element_text(size = 15, hjust = 0.5)) 
x
data.survdiff <- survdiff(Surv(OS.time,OS)~ group_B ,data = LUAD_TPM_clin)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
ci <- paste0(sprintf("%.3f",HR)," [",sprintf("%.3f",low95),", ",sprintf("%.3f",up95),"]")

x$plot <- x$plot + 
  annotate("text", x = 1470, y = 0.83, 
           label = paste0("\n HR (95% CI) = ",ci),  
           size = 5, color = "black", hjust = 0)+
  theme(text = element_text(size = 15))
x$plot


########################################################################################################

#LM high vs LM low DESeq2-DEG
library(DESeq2)
library(edgeR)
library(limma)
library(pheatmap)
library(BiocParallel)
library(RColorBrewer)
library(IHW)
library(vsn)
library(ggplot2)
library(genefilter)
library(SummarizedExperiment)
library(VennDiagram)
library(vsn)
library(reshape2)
library(apeglm)
library(clusterProfiler)
library(org.Hs.eg.db)


setwd("G:/LUAD_scRNASeq/LUAD_bulk_RNASeq/")
LUAD_count = read.csv("G:/MR_mGWAS/LUAD_RNASeq/TCGA_LUAD_count_symbol.csv",row.names = 1,check.names = F)
LUAD_count = LUAD_count[,order(colnames(LUAD_count))]
colnames(LUAD_count)=substr(colnames(LUAD_count),1,15)
sum(duplicated(colnames(LUAD_count)))
LUAD_count = LUAD_count[,!duplicated(colnames(LUAD_count))]
LUAD_count = LUAD_count[rownames(LUAD_TPM),LUAD_TPM_clin$sample]
colnames(LUAD_count)=substr(colnames(LUAD_count),1,12)

myfactor=LM_gsva
myfactor$sample = rownames(myfactor)
myfactor$condition = myfactor$LM_group
myfactor=myfactor[,-c(1,2)]
myfactor$condition = factor(myfactor$condition,levels = c("high_LM","low_LM"))
myfactor = myfactor[order(myfactor$condition),]
LUAD_count = LUAD_count[,rownames(myfactor)]

count_matrix = as.matrix(LUAD_count)
count_Matrix=count_matrix
register(SnowParam(workers=4, type = "SOCK")) 
count_DESeq<-DESeqDataSetFromMatrix(countData = count_Matrix,colData = myfactor,design=~condition)

count_DESeq <- DESeq(count_DESeq)
print(resultsNames(count_DESeq))

res<-results(count_DESeq,contrast = c("condition","high_LM","low_LM"))
print(res)
print(summary(res)) 
res=as.data.frame(res)
res=res[order(res$padj),]
res=na.omit(res)
write.csv(res,"G:/LUAD_scRNASeq/LUAD_bulk_RNASeq/LUAD_LM_DEG",row.names = T)

#vp
res = read.csv("G:/LUAD_scRNASeq/LUAD_bulk_RNASeq/LUAD_LM_DEG",header = T,row.names = 1)
library(ggplot2)
library(ggrepel)  
library(pheatmap)
library(export)
data = res
data$Name = rownames(data)
data <- na.omit(data)

log2FC = 1 
pvalue = 0.01

data$sig[(-1*log10(data$padj) < -1*log10(pvalue)|data$padj=="NA")| data$log2FoldChange < log2FC & data$log2FoldChange > -(log2FC)] <- "NotSig"
data$sig[-1*log10(data$padj) >= -1*log10(pvalue) & data$log2FoldChange >= log2FC] <- "Up"
data$sig[-1*log10(data$padj) >= -1*log10(pvalue) & data$log2FoldChange <= -(log2FC)] <- "Down"

table(data$sig)
this_tile <- paste0('Cutoff for padj is ',pvalue,
                    '\nThe number of up gene is ',nrow(data[data$sig =="Up",]) ,
                    '\nThe number of down gene is ',nrow(data[data$sig =="Down",]))                 
this_tile

PvalueLimit = 1.792274e-44
data$label=ifelse(data$padj <= PvalueLimit & data$sig != "NotSig" , as.character(data$Name), '')
write.csv(data,"G:/LUAD_scRNASeq/Fig_and_chart/Fig2C_DEG_LM_high_vs_LM_low.csv")

# 绘图
x = ggplot(data,aes(log2FoldChange,-1*log10(padj))) +    
  geom_point(aes(color = sig)) +                          
  xlab("log2 fold change") + 
  ylab("-log10 p-value") +
  ggtitle( this_tile ) +                             
  # scale_color_manual(values = c("blue","green","red")) + 
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  geom_hline(yintercept=-log10(pvalue),linetype=2)+        
  geom_vline(xintercept=c(-(log2FC),log2FC),linetype=2)+ 
  geom_text_repel(aes(x = log2FoldChange,                   
                      y = -1*log10(padj),          
                      label=label),                       
                  max.overlaps = 10000,                    
                  size=3.2,                                  
                  box.padding=unit(0.5,'lines'),          
                  point.padding=unit(0.1, 'lines'), 
                  segment.color='black',      
                  show.legend=FALSE)                     

x

#################################################################################################

#ESTIMATE,CIBERSORT,CIBERSORTx(online)
setwd("G:/MR_mGWAS/LUAD_RNASeq/")
LUAD_group = LUAD_TPM_clin[,c("X_PATIENT","risk_B","group_B")]
colnames(LUAD_group) = c("Row.names","risk","group")
LUAD_group$LM_group = LM_gsva$LM_group[match(LUAD_group$Row.names,rownames(LM_gsva))]

library(estimate)
library(dplyr)

LUAD_TPM = t(LUAD_TPM_clin[,12:30492])
colnames(LUAD_TPM) = LUAD_TPM_clin$X_PATIENT
LUAD_TPM = LUAD_TPM[,colnames(LUAD_TPM) %in% LUAD_group$Row.names]
LUAD_TPM = LUAD_TPM[,order(colnames(LUAD_TPM))]
all(colnames(LUAD_TPM) == LUAD_group$Row.names)
write.table(LUAD_TPM,"LUAD_TPM.txt",sep="\t",row.names = T,quote = F)
write.csv(LUAD_group,"LUAD_risk_group.csv",row.names = F)
gc()

estimate <- function(pro){
  input.f=paste0(pro,'_TPM.txt')
  output.f=paste0(pro,'_estimate_gene.gct')
  output.ds=paste0(pro,'_estimate_score.gct')
  library(estimate)
  filterCommonGenes(input.f=input.f,
                    output.f=output.f ,
                    id="GeneSymbol")
  estimateScore(input.ds = output.f,
                output.ds=output.ds,
                platform="illumina")   ## 注意platform,RNA-Seq为illumina，芯片为affy
  scores=read.table(output.ds,skip = 2,header = T)
  rownames(scores)=scores[,1]
  scores=t(scores[,3:ncol(scores)])
  return(scores)
}
pro='LUAD'
ESTIMATE_scores=estimate(pro)


LUAD_TPM = read.table("LUAD_TPM.txt",sep="\t",header = T)
LUAD_group = read.csv("LUAD_risk_group.csv",header = T)
library(e1071)
#devtools::install_github("Moonerss/CIBERSORT")
library(CIBERSORT)
data(LM22)
LUAD_TPM_mat = as.matrix(LUAD_TPM)
CIBER_res = cibersort(sig_matrix = LM22, mixture_file =LUAD_TPM_mat)


#boxplot，violin plot
rownames(ESTIMATE_scores) = gsub("\\.", "-", rownames(ESTIMATE_scores))
all(rownames(ESTIMATE_scores)==LUAD_group$Row.names)
ESTIMATE_scores = data.frame(ESTIMATE_scores)
ESTIMATE_scores$Group = LUAD_group$LM_group
library(tidyr)
df_long <- pivot_longer(data.frame(ESTIMATE_scores), cols = c(StromalScore,ImmuneScore,ESTIMATEScore), names_to = "Variable", values_to = "Value")
library(ggplot2)
library(ggpubr)

p <- ggplot(df_long, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +  
  facet_wrap(~Variable,scales = "free") +
  stat_compare_means(aes(label = ..p.signif..), method = "wilcox.test",hjust = -0.8)+
  scale_fill_brewer(palette = "Pastel1") +  
  theme_minimal() +  
  labs(title = "ESTIMATE", x = NULL, y = "Scores")  
print(p)


#CIBERSORT
rownames(CIBER_res) = gsub("\\.", "-", rownames(CIBER_res))
all(rownames(CIBER_res)==LUAD_group$Row.names)
CIBER_res = data.frame(CIBER_res)
CIBER_res$Group = LUAD_group$LM_group
CIBER_res = CIBER_res[,!colnames(CIBER_res) %in% c("P.value","Correlation","RMSE")]
library(tidyr)
df_long <- pivot_longer(data.frame(CIBER_res), cols = -ncol(data.frame(CIBER_res)), names_to = "Variable", values_to = "Value")
library(ggplot2)
library(ggpubr)
p <- ggplot(df_long, aes(x = Variable, y = Value, fill = Group)) +
  geom_boxplot(outlier.shape = NA) + 
  scale_fill_brewer(palette = "Pastel1") + 
  theme_minimal() +  
  labs(title = "CIBERSORT", x = "immune cell type", y = "Porp")+  
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
p <- p + stat_compare_means(aes(group = Group), method = "t.test", label = "p.signif",hide.ns = T)
print(p)

median_score = ESTIMATE_scores %>% group_by(Group) %>%
  summarise(StromalScore_median = median(StromalScore),ImmuneScore_median = median(ImmuneScore),ESTIMATEScore_median = median(ESTIMATEScore))
p_EST = data.frame(Score = c("StromalScore","ImmuneScore","ESTIMATEScore"),
                   pvalue = c(
                     wilcox.test(StromalScore ~ Group,data = ESTIMATE_scores)$p.value,
                     wilcox.test(ImmuneScore ~ Group,data = ESTIMATE_scores)$p.value,
                     wilcox.test(ESTIMATEScore ~ Group,data = ESTIMATE_scores)$p.value
                   ))

median_CIB = CIBER_res %>% group_by(Group) %>% summarise(across(1:22,median))
p_CIB = map_dfr(1:22,function(i){
  p_value = wilcox.test(CIBER_res[[i]] ~ CIBER_res$Group)$p.value
  data.frame(Score = colnames(CIBER_res)[i],p_value=p_value)
})


############################################################################
##univ/Multi-Cox

library(GSVA)
library(survival)
library(DT)
library(survminer)
library(stringr)
library(dplyr)
setwd("G:/LUAD_scRNASeq/")
pheo = read.table("G:/LUAD_scRNASeq/LUAD_pheo_survival/CDR_LUAD_patient_pheo.txt",sep="\t",check.names = F,header = T)

LUAD_TPM_clin = read.csv("G:/MR_mGWAS/LUAD_RNASeq/LUAD_TPM_Clin.csv",header = T,
                         row.names = 1,check.names = F)
LUAD_count = read.csv("G:/MR_mGWAS/LUAD_RNASeq/TCGA_LUAD_count_symbol.csv",row.names = 1,check.names = F)
LUAD_TPM = t(LUAD_TPM_clin[,12:ncol(LUAD_TPM_clin)])
colnames(LUAD_TPM) = LUAD_TPM_clin$X_PATIENT
load("cc_Wu_Zhou_2021_LMn.RData")
LMn = intersect(LMn$gene,rownames(LUAD_TPM))
LMn = list(LM = LMn)

#LMn GSVA
LM_gsva <- gsva(expr = LUAD_TPM, gset.idx.list = LMn, kcdf = "Gaussian",parallel.sz = 8)
LM_gsva = as.data.frame(t(LM_gsva))
LM_gsva$LM_group = ifelse(LM_gsva$LM > median(LM_gsva$LM),"high_LM","low_LM")
all(rownames(LM_gsva) == LUAD_TPM_clin$X_PATIENT)
LUAD_TPM_clin$LM_group = LM_gsva$LM_group

LUAD_TPM_clin = LUAD_TPM_clin[!(LUAD_TPM_clin$OS ==0 & LUAD_TPM_clin$OS.time<30),]
LUAD_TPM_clin = LUAD_TPM_clin[!is.na(LUAD_TPM_clin$OS.time),]

library(timeROC)
fit_B <- coxph(Surv(OS.time,OS)~ CCT6A+KRT8+LYPD3+HMGA1+PSMB7 ,data = LUAD_TPM_clin)
risk_B = predict(fit_B,newx=LUAD_TPM_clin[,c("CCT6A","KRT8","LYPD3","HMGA1","PSMB7")],type="lp")

LUAD_TPM_clin$risk_B = risk_B
LUAD_TPM_clin$group_B = ifelse(LUAD_TPM_clin$risk_B>median(LUAD_TPM_clin$risk_B),"high-risk","low-risk")
surv_fit_B = survfit(Surv(OS.time,OS)~ group_B ,data = LUAD_TPM_clin)

cols = c("sampleID","additional_pharmaceutical_therapy","additional_radiation_therapy","additional_surgery_locoregional_procedure",
         "additional_surgery_metastatic_procedure","age_at_initial_pathologic_diagnosis","histological_type",
         "history_of_neoadjuvant_treatment","pathologic_M","pathologic_N","pathologic_T","pathologic_stage",
         "tobacco_smoking_history","tobacco_smoking_history_indicator")
pheo = pheo[pheo$sampleID %in% LUAD_TPM_clin$sample,cols]

cols = c("sampleID","age_at_initial_pathologic_diagnosis","pathologic_M","pathologic_N","pathologic_T","pathologic_stage")
pheo = pheo[,cols]
colnames(pheo)=c("sampleID","age","M_stage","N_stage","T_stage","stage")
pheo = pheo[pheo$N_stage != "" & pheo$T_stage != "" & pheo$M_stage != "" & pheo$stage !="",]
pheo = pheo[pheo$stage != "[Discrepancy]",]

LUAD_TPM_clin = LUAD_TPM_clin[LUAD_TPM_clin$sample %in% pheo$sampleID,]
LUAD_TPM_clin = LUAD_TPM_clin[order(LUAD_TPM_clin$sample),]
pheo = pheo[order(pheo$sampleID),] 

LUAD_TPM_clin = cbind(LUAD_TPM_clin,pheo)
clin = LUAD_TPM_clin

clin$age_group = ifelse(clin$age>60,">60","<60")
clin$M_group = ifelse(clin$M_stage != "M0",">M0","M0")
clin$T_group = ifelse(clin$T_stage > "T2",">T2","<=T2")
clin$N_group = ifelse(clin$N_stage > "N0",">N0","N0")
clin$stage_group = ifelse(clin$stage > "Stage II",">II","<=II")
clin$risk_group = ifelse(clin$risk_B > median(clin$risk_B),"high-risk","low-risk")

clin$age_group = factor(clin$age_group, levels = c("<60", ">60"))
clin$M_group = factor(clin$M_group, levels = c("M0", ">M0"))
clin$T_group = factor(clin$T_group, levels = c("<=T2", ">T2"))
clin$N_group = factor(clin$N_group, levels = c("N0", ">N0"))
clin$stage_group = factor(clin$stage_group, levels = c("<=II", ">II"))
clin$risk_group <- factor(clin$group, levels = c("low-risk", "high-risk"))

fit1=coxph(Surv(OS.time,OS)~ age_group ,data = clin) 
fit2=coxph(Surv(OS.time,OS)~ M_group ,data = clin)
fit3=coxph(Surv(OS.time,OS)~ T_group ,data = clin)
fit4=coxph(Surv(OS.time,OS)~ N_group ,data = clin)
fit5=coxph(Surv(OS.time,OS)~ stage_group ,data = clin)
fit6=coxph(Surv(OS.time,OS)~ risk_group ,data = clin)

unicox=data.frame(matrix(0,nrow=6,ncol=5))
rownames(unicox) = c("age","M_stage","T_stage","N_stage","stage","LM_risk")
colnames(unicox) = c("","HR","lower.95","upper.95","pvalue")
unicox[,1] = c(">60/<60",">M0/M0",">T2/<=T2",">N0/N0",">II/<=II","High_risk/Low_risk")

for(i in 1:6){
  a = get(paste0("fit",as.character(i)))
  summ = summary(a)
  unicox[i,2] = summ$coefficients[,"exp(coef)"]
  unicox[i,3] = summ$conf.int[,"lower .95"]
  unicox[i,4] = summ$conf.int[,"upper .95"]
  unicox[i,5] = summ$coefficients[,"Pr(>|z|)"]
}

fit7=coxph(Surv(OS.time,OS)~ age_group+M_group+T_group+N_group+stage_group+risk_group ,data = clin)
summary(fit7)
Mulcox=data.frame(matrix(0,nrow=6,ncol=5))
rownames(Mulcox) = c("age","M_stage","T_stage","N_stage","stage","LM_risk")
colnames(Mulcox) = c("","HR","lower.95","upper.95","pvalue")
Mulcox[,1] = c(">60/<60",">M0/M0",">T2/<=T2",">N0/N0",">II/<=II","High_risk/Low_risk")
Mulcox[,2] = summary(fit7)$coefficients[,"exp(coef)"]
Mulcox[,3] = summary(fit7)$conf.int[,"lower .95"]
Mulcox[,4] = summary(fit7)$conf.int[,"upper .95"]
Mulcox[,5] = summary(fit7)$coefficients[,"Pr(>|z|)"]

unicox = cbind(rownames(unicox),unicox)
Mulcox = cbind(rownames(Mulcox),Mulcox)
colnames(unicox)[1] = "Variable"
colnames(Mulcox)[1] = "Variable"


Mulcox[,1:5] <- data.frame(lapply(Mulcox[,1:5], function(x) if(is.numeric(x)) round(x, 2) else x))
unicox[,1:5] <- data.frame(lapply(unicox[,1:5], function(x) if(is.numeric(x)) round(x, 2) else x))
unicox <- unicox %>%
  mutate(pvalue = ifelse(round(pvalue, 3) < 0.001, "<0.001", sprintf("%.3f", round(pvalue, 3))))
Mulcox <- Mulcox %>%
  mutate(pvalue = ifelse(round(pvalue, 3) < 0.001, "<0.001", sprintf("%.3f", round(pvalue, 3))))

unicox <- unicox %>%
  mutate(formatted_column = paste0(HR, " (", lower.95, " - ", upper.95, ")"))
unicox$HR = unicox$formatted_column
unicox = unicox[,c("Variable","Var.2","HR","pvalue")]
Mulcox <- Mulcox %>%
  mutate(formatted_column = paste0(HR, " (", lower.95, " - ", upper.95, ")"))
Mulcox$HR = Mulcox$formatted_column
Mulcox = Mulcox[,c("Variable","Var.2","HR","pvalue")]


