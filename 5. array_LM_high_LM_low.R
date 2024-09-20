library(GSVA)
library(survival)
library(DT)
library(survminer)
library(stringr)
library(dplyr)

LM_B = c("CCT6A", "KRT8", "LYPD3", "HMGA1", "PSMB7")
setwd("E:/Bioinfor_reproduction/Rloop/microarray/LUAD")
load("./GSE26939/exprSet.Rdata")
load("./GSE26939/pd.Rdata")

all(LM_B %in% rownames(exprSet))
exprSet = exprSet[,rownames(pd)[pd$Histology=="TP"]]
pd = pd[pd$Histology=="TP",]

clin_exp = data.frame(cbind(pd,t(exprSet)))
library(ROCR)
library(caret)
library(ggplot2)
library(survminer)
library(car)

fit_B <- coxph(Surv(OS.Time,OS)~ CCT6A+KRT8+LYPD3+HMGA1+PSMB7 ,data = clin_exp)
library(timeROC)
risk_B = predict(fit_B,newx=clin_exp[,c("CCT6A","KRT8","LYPD3","HMGA1","PSMB7")],type="lp")

clin_exp$risk_B = risk_B
clin_exp$group_B = ifelse(clin_exp$risk_B>median(clin_exp$risk_B),"high-risk","low-risk")

surv_fit_B = survfit(Surv(OS.Time,OS)~ group_B ,data = clin_exp)
x=ggsurvplot(surv_fit_B,
             pval = TRUE,
             main = "Survival curve",
             conf.int = F,
             break.time.by = 1,
             pval.coord = c(4, 0.9), 
             xlim = c(0, 10),
             risk.table = TRUE, # Add risk table
             risk.table.col = "strata", # Change risk table color by groups
             tables.height = 0.2, 
             tables.theme = theme_cleantable(), 
             linetype = "strata", # Change line type by groups
             surv.median.line = "hv" # Specify median survival
)
x$plot = x$plot + ggtitle("5 LM risk gene in GSE26939") + theme(plot.title = element_text(size = 15, hjust = 0.5)) 
x

data.survdiff <- survdiff(Surv(OS.Time,OS)~ group_B ,data = clin_exp)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
ci <- paste0(sprintf("%.3f",HR)," [",sprintf("%.3f",low95),", ",sprintf("%.3f",up95),"]")

x$plot <- x$plot + 
  annotate("text", x = 3.8, y = 0.83, 
           label = paste0("\n HR (95% CI) = ",ci),  
           size = 5, color = "black", hjust = 0)+
  theme(text = element_text(size = 15))

setwd("E:/Bioinfor_reproduction/Rloop/microarray/LUAD")
load("./GSE30219/exprSet.Rdata")
load("./GSE30219/pd.Rdata")

all(LM_B %in% rownames(exprSet))
#去除NT
exprSet = exprSet[,rownames(pd)[pd$Histology=="TP"]]
pd = pd[pd$Histology=="TP",]

clin_exp = data.frame(cbind(pd,t(exprSet)))
library(ROCR)
library(caret)
library(ggplot2)
library(survminer)
library(car)

fit_B <- coxph(Surv(OS.Time,OS)~ CCT6A+KRT8+LYPD3+HMGA1+PSMB7 ,data = clin_exp)
library(timeROC)
risk_B = predict(fit_B,newx=clin_exp[,c("CCT6A","KRT8","LYPD3","HMGA1","PSMB7")],type="lp")

clin_exp$risk_B = risk_B
clin_exp$group_B = ifelse(clin_exp$risk_B>median(clin_exp$risk_B),"high-risk","low-risk")

surv_fit_B = survfit(Surv(OS.Time,OS)~ group_B ,data = clin_exp)
x=ggsurvplot(surv_fit_B,
             pval = TRUE,
             main = "Survival curve",
             conf.int = F,
             break.time.by = 1,
             pval.coord = c(4, 0.3), 
             xlim = c(0, 10),
             risk.table = TRUE, # Add risk table
             risk.table.col = "strata", # Change risk table color by groups
             tables.height = 0.2, 
             tables.theme = theme_cleantable(), 
             linetype = "strata", # Change line type by groups
             surv.median.line = "hv" # Specify median survival
)
x$plot = x$plot + ggtitle("LM risk gene in GSE30129") + theme(plot.title = element_text(size = 15, hjust = 0.5)) 
x

data.survdiff <- survdiff(Surv(OS.Time,OS)~ group_B ,data = clin_exp)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
ci <- paste0(sprintf("%.3f",HR)," [",sprintf("%.3f",low95),", ",sprintf("%.3f",up95),"]")

x$plot <- x$plot + 
  annotate("text", x = 3.8, y = 0.23, 
           label = paste0("\n HR (95% CI) = ",ci),   ###添加P和HR 95%CI
           size = 5, color = "black", hjust = 0)+
  theme(text = element_text(size = 15))
x$plot

load("G:/LUAD_scRNASeq/cc_Wu_Zhou_2021_LMn.RData")
LMn = intersect(LMn$gene,rownames(exprSet))
LMn = list(LM = LMn)

load("./GSE13213/exprSet.Rdata")
load("./GSE13213/pd.Rdata")

all(LM_B %in% rownames(exprSet))
table(pd$Histology)
all(colnames(exprSet) == pd$Sample)
clin_exp = data.frame(cbind(pd,t(exprSet)))

LM_gsva <- gsva(expr = as.matrix(exprSet), gset.idx.list = LMn, kcdf = "Gaussian",parallel.sz = 8)
LM_gsva = as.data.frame(t(LM_gsva))
LM_gsva$LM_group = ifelse(LM_gsva$LM > median(LM_gsva$LM),"high_LM","low_LM")
all(rownames(LM_gsva) == clin_exp$Sample)
clin_exp$LM_group = LM_gsva$LM_group

fit <- survfit(Surv(OS.Time, OS) ~ LM_group, data = clin_exp)
x=ggsurvplot(fit,
             pval = TRUE,
             pval.coord = c(0.5, 0.25),
             main = "Survival curve",
             conf.int = F,
             break.time.by = 1,
             xlim = c(0, 10),
             risk.table = TRUE, # Add risk table
             risk.table.col = "strata", # Change risk table color by groups
             tables.height = 0.2,
             tables.theme = theme_cleantable(),
             linetype = "strata", # Change line type by groups
             surv.median.line = "hv" # Specify median survival
)
x$plot = x$plot + ggtitle("GSE13123") + theme(plot.title = element_text(size = 15, hjust = 0.5)) 
x

data.survdiff <- survdiff(Surv(OS.Time,OS)~ LM_group ,data = clin_exp)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
ci <- paste0(sprintf("%.3f",HR)," [",sprintf("%.3f",low95),", ",sprintf("%.3f",up95),"]")

x$plot <- x$plot + 
  annotate("text", x = 0.3, y = 0.18, 
           label = paste0("\n HR (95% CI) = ",ci),  
           size = 5, color = "black", hjust = 0)+
  theme(text = element_text(size = 15))
x$plot


load("./GSE31210/exprSet.Rdata")
load("./GSE31210/pd.Rdata")

all(LM_B %in% rownames(exprSet))
table(pd$Histology)
exprSet = exprSet[,rownames(pd)[pd$Histology=="TP"]]
pd = pd[pd$Histology=="TP",]
all(colnames(exprSet) == pd$Sample)
clin_exp = data.frame(cbind(pd,t(exprSet)))

LM_gsva <- gsva(expr = as.matrix(exprSet), gset.idx.list = LMn, kcdf = "Gaussian",parallel.sz = 8)
LM_gsva = as.data.frame(t(LM_gsva))
colnames(LM_gsva) = "LM"
LM_gsva$LM_group = ifelse(LM_gsva$LM > median(LM_gsva$LM),"high_LM","low_LM")
all(rownames(LM_gsva) == clin_exp$Sample)
clin_exp$LM_group = LM_gsva$LM_group


fit <- survfit(Surv(OS.Time, OS) ~ LM_group, data = clin_exp)
x=ggsurvplot(fit,
             pval = TRUE,
             pval.coord = c(0.5, 0.25),
             main = "Survival curve",
             conf.int = F,
             break.time.by = 1,
             xlim = c(0, 10),
             risk.table = TRUE, # Add risk table
             risk.table.col = "strata", # Change risk table color by groups
             tables.height = 0.2, 
             tables.theme = theme_cleantable(), 
             linetype = "strata", # Change line type by groups
             surv.median.line = "hv" # Specify median survival
)
x$plot = x$plot + ggtitle("GSE31210") + theme(plot.title = element_text(size = 15, hjust = 0.5)) 
x

data.survdiff <- survdiff(Surv(OS.Time,OS)~ LM_group ,data = clin_exp)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
ci <- paste0(sprintf("%.3f",HR)," [",sprintf("%.3f",low95),", ",sprintf("%.3f",up95),"]")

x$plot <- x$plot + 
  annotate("text", x = 0.3, y = 0.18, 
           label = paste0("\n HR (95% CI) = ",ci),   
           size = 5, color = "black", hjust = 0)+
  theme(text = element_text(size = 15))
x$plot



