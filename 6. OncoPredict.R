setwd("G:/LUAD_scRNASeq/LUAD_bulk_RNASeq/")

#LUAD_TPM是log2(TPM+1)
#LUAD_TPM = 2^LUAD_TPM-1
LUAD_TPM_clin = read.csv("G:/MR_mGWAS/LUAD_RNASeq/LUAD_TPM_Clin.csv",header = T,
                         row.names = 1,check.names = F)
LUAD_count = read.csv("G:/MR_mGWAS/LUAD_RNASeq/TCGA_LUAD_count_symbol.csv",row.names = 1,check.names = F)
LUAD_TPM = t(LUAD_TPM_clin[,12:ncol(LUAD_TPM_clin)])
dim(LUAD_TPM)

library(oncoPredict)
library(reshape2)
library(ggpubr)
#读取CTRP和GDSC的训练数据
dir='G:/MR_mGWAS/LUAD_RNASeq/DrugResponse_data/DataFiles/DataFiles/Training Data/'
th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
GDSC2_Expr = readRDS(file=file.path(dir,'GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
dim(GDSC2_Expr)  
GDSC2_Expr[1:4, 1:4]
boxplot(GDSC2_Expr[,1:4])
df=melt(GDSC2_Expr[,1:4])
head(df)
p1=ggboxplot(df, "Var2", "value") +th

# Read GDSC2 response data. rownames() are samples, colnames() are drugs. 
dir
GDSC2_Res = readRDS(file = file.path(dir,"GDSC2_Res.rds"))
dim(GDSC2_Res)  # 805 198
GDSC2_Res[1:4, 1:4]
p2=ggboxplot(melt(GDSC2_Res[,1:4]), "Var2", "value") +th ; p2
# IMPORTANT note: here I do e^IC50 since the IC50s are actual ln values/log transformed already, and the calcPhenotype function Paul #has will do a power transformation (I assumed it would be better to not have both transformations)
GDSC2_Res <- exp(GDSC2_Res)  
p3=ggboxplot(melt(GDSC2_Res[,1:4]), "Var2", "value") +th ; p3

library(patchwork)
p1+p2+p3

#不同cell-line IC50最高和最低的药物
IC50_sort_GDSC = round(apply(GDSC2_Res, 1, function(x){
  return(c(
    head(sort(x)),
    tail(sort(x))
  ))
}),2)

drug_sort_GDSC = apply(GDSC2_Res, 1, function(x){ 
  names(x)=gsub('_[0-9]*','',colnames(GDSC2_Res))
  return(c(
    names(head(sort(x))),
    names(tail(sort(x)))
  ))
})

#GDSC的表达矩阵是log转化的芯片表达，CTRP2为没有对数的TPM/CPM
#GDSC1: 367 drugs,958 cell-types
GDSC1_Expr = readRDS(file=file.path(dir,'GDSC1_Expr (RMA Normalized and Log Transformed).rds'))
GDSC1_Res = readRDS(file = file.path(dir,"GDSC1_Res.rds"))
GDSC1_Res <- exp(GDSC1_Res)  
#CTRP:829 drugs, 545 cell-types
CTRP2_Expr = readRDS(file=file.path(dir,'CTRP2_Expr (TPM, not log transformed).rds'))
CTRP2_Res = readRDS(file = file.path(dir,"CTRP2_Res.rds"))

#predict IC50
CTRP2_LUAD = calcPhenotype(trainingExprData = CTRP2_Expr,
                           trainingPtype = CTRP2_Res,
                           testExprData = as.matrix(2^LUAD_TPM-1),#需要matrix,和CTRP2训练数据量纲一致
                           batchCorrect = 'eb',  
                           powerTransformPhenotype = F,
                           minNumSamples = 20,
                           printOutput = T,
                           removeLowVaryingGenes = 0.3,
                           removeLowVaringGenesFrom = "homogenizeData",
                           folder=F)
GDSC1_LUAD = calcPhenotype(trainingExprData = GDSC1_Expr,
                           trainingPtype = GDSC1_Res,
                           testExprData = as.matrix(LUAD_TPM),#需要matrix,log2TPM
                           batchCorrect = 'eb',  
                           powerTransformPhenotype = T,
                           minNumSamples = 20,
                           printOutput = T,
                           removeLowVaryingGenes = 0.2,
                           removeLowVaringGenesFrom = "homogenizeData",
                           folder=F)
GDSC2_LUAD = calcPhenotype(trainingExprData = GDSC2_Expr,
                           trainingPtype = GDSC2_Res,
                           testExprData = as.matrix(LUAD_TPM),#需要matrix,log2TPM
                           batchCorrect = 'eb',  
                           powerTransformPhenotype = T,
                           minNumSamples = 20,
                           printOutput = T,
                           removeLowVaryingGenes = 0.2,
                           removeLowVaringGenesFrom = "homogenizeData",
                           folder=F)
setwd("G:/LUAD_scRNASeq/LUAD_bulk_RNASeq/")
write.csv(GDSC1_LUAD,"GDSC1_LUAD_IC50_predict.csv",row.names = T)
write.csv(GDSC2_LUAD,"GDSC2_LUAD_IC50_predict.csv",row.names = T)
write.csv(CTRP2_LUAD,"CTRP2_LUAD_IC50_predict.csv",row.names = T)

#risk score和IC50的相关系数
library(GSVA)
GDSC1_LUAD = read.csv("GDSC1_LUAD_IC50_predict.csv",header=T)
load("G:/LUAD_scRNASeq/cc_Wu_Zhou_2021_LMn.RData")
LMn = intersect(LMn$gene,rownames(LUAD_TPM))
LMn = list(LM = LMn)
LM_gsva <- gsva(expr = LUAD_TPM, gset.idx.list = LMn, kcdf = "Gaussian",parallel.sz = 8)
LM_gsva = as.data.frame(t(LM_gsva))
LM_gsva$LM_group = ifelse(LM_gsva$LM > median(LM_gsva$LM),"high_LM","low_LM")
LUAD_group = LM_gsva

res_GDSC1 = data.frame(matrix(0,nrow = 2,ncol=ncol(GDSC1_LUAD)))
rownames(res_GDSC1) = c("cor","p")
colnames(res_GDSC1)=colnames(GDSC1_LUAD)
all(rownames(GDSC1_LUAD) == LUAD_group$Row.names)
for(i in 1:ncol(res_GDSC1)){
  cor_test = cor.test(LUAD_group$LM,GDSC1_LUAD[,i])
  res_GDSC1[1,i]=round(cor_test$estimate,3)
  res_GDSC1[2,i]=cor_test[["p.value"]]
}
res_GDSC1 = data.frame(t(res_GDSC1))
res_GDSC1 = res_GDSC1[order(res_GDSC1$p,decreasing = F),]

#risk和IC50相关性散点图
drugs = c(rownames(res_GDSC1)[res_GDSC1$cor<0][1:4],rownames(res_GDSC1)[res_GDSC1$cor>0][1:4])
df = cbind(LUAD_group,GDSC1_LUAD)

res = res_GDSC1[rownames(res_GDSC1) %in% drugs,]
write.csv(res,"G:/LUAD_scRNASeq/Fig_and_chart/Fig2E.csv",row.names = T)

pic_list = list()
for(i in 1:8){
  drug = drugs[i]
  cor_value = res_GDSC1[drug,1]
  p_value = res_GDSC1[drug,2]
  cor_label = paste0("r=",cor_value,"\np = ",format(p_value,digits=3,scientific = T))
  
  p = ggplot(df,aes(x = LM,y=!!sym(drug))) + 
    geom_point(color = "blue",alpha = 0.6) + 
    geom_smooth(method = "lm",color = "red",se=FALSE) + 
    theme_minimal() + 
    annotate("text",x=Inf,y=Inf,label = cor_label,hjust = 1.1,vjust = 1.1,size =2,color = "black") + 
    labs(x="LM score",y=drug) +
    theme(plot.title = element_text(size=10),axis.title.x = element_text(size=8),axis.title.y = element_text(size=8))
  pic_list[[i]] = p
}

combined_plot <- ggarrange(
  pic_list[[1]], pic_list[[2]], pic_list[[3]], pic_list[[4]], pic_list[[5]], pic_list[[6]], pic_list[[7]], pic_list[[8]], 
  ncol = 4, nrow = 2,
  common.legend = TRUE, legend = "right"
)

ann_plot = annotate_figure(combined_plot,top = text_grob("cor of IC50 vs LM score",size = 14,face="bold"))
print(ann_plot)

