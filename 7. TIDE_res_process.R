
library(RColorBrewer)

setwd("G:/LUAD_scRNASeq/LUAD_bulk_RNASeq/")
LUAD_TPM = read.table("LUAD_TPM.txt",sep="\t",header = T)
load("G:/LUAD_scRNASeq/cc_Wu_Zhou_2021_LMn.RData")
LMn = intersect(LMn$gene,rownames(LUAD_TPM))
LMn = list(LM = LMn)
LUAD_TPM = as.matrix(LUAD_TPM)
LM_gsva <- gsva(expr = LUAD_TPM, gset.idx.list = LMn, kcdf = "Gaussian",parallel.sz = 8)
LM_gsva = as.data.frame(t(LM_gsva))
LM_gsva$LM_group = ifelse(LM_gsva$LM > median(LM_gsva$LM),"high_LM","low_LM")

TIDE_expr = t(apply(LUAD_TPM,1,function(x) {x-mean(x)}))
write.table(TIDE_expr,"TIDE.txt",sep = "\t",row.names = T,quote = F)

LUAD_group = LM_gsva
TIDE_res_NSCLC = read.csv("TIDE_imm_res_NSCLC.csv",header = T,check.names = F)
TIDE_res_Mela = read.csv("TIDE_imm_res_melanoma.csv",header = T,check.names = F)
rownames(TIDE_res_NSCLC) = TIDE_res_NSCLC$Patient
TIDE_res_NSCLC = TIDE_res_NSCLC[rownames(LUAD_group),]
LUAD_group = cbind(LUAD_group,TIDE_res_NSCLC[,2:ncol(TIDE_res_NSCLC)])


col_type = sapply(LUAD_group[,4:ncol(LUAD_group)],is.numeric)
pvalue_res = data.frame(matrix(0,nrow=1,ncol=ncol(LUAD_group)-3))
colnames(pvalue_res) = colnames(LUAD_group)[4:ncol(LUAD_group)]
for(i in 1:(ncol(LUAD_group)-3)){
  df = list()
  if(col_type[i]){
    df$group = LUAD_group$LM_group
    df$value = LUAD_group[,i+3]
    test = wilcox.test(value ~ group,data=df)
    pvalue_res[1,i]=test$p.value
  }else{
    df$group = LUAD_group[,i+3]
    df$value = LUAD_group$LM
    test = wilcox.test(value ~ group,data=df)
    pvalue_res[1,i]=test$p.value
  }
}

library(ggplot2)
library(RColorBrewer)
library(ggsignif)
library(ggpubr)

#responser and No.benefit vs LM score
p = ggplot(LUAD_group,aes(x = Responder,y=LM,fill = Responder)) + 
  geom_boxplot(width =0.35,position = position_dodge(width = 0.2)) + 
  scale_fill_manual(values = brewer.pal(n=3,name = "Set2")) + 
  theme_classic() + 
  theme(aspect.ratio = 0.9,axis.text = element_text(size = 12),axis.title = element_text(size = 14),legend.position = "none") +
  labs(x="TIDE responder group",y="LM score",title = "       LM score vs responder group") + 
  geom_signif(comparisons = list(c("True","False")),map_signif_level = T,
              y_position = max(LUAD_group$LM)+0.2,annotations = "*")
p


library(ggplot2)
library(ggpubr)

df = LUAD_group[,c("LM_group","IFNG","CD8","Exclusion","TAM M2")]
df_long <- tidyr::pivot_longer(df, cols = c("IFNG","CD8","Exclusion","TAM M2"), 
                               names_to = "Variable", values_to = "Value")

plot_box <- function(data, x, y, group_var) {
  ggboxplot(data, x = group_var, y = y, color = group_var, palette = "jco") +
    stat_compare_means(aes(label = ..p.signif..), method = "t.test",label.x=1.4) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45, vjust = 0.7))
}

plots <- list()
for (variable in unique(df_long$Variable)) {
  plots[[variable]] <- plot_box(df_long[df_long$Variable == variable, ], "LM_group", "Value", "LM_group") +
    ggtitle(variable)
}

combined_plot <- ggarrange(
  plots[[1]], plots[[2]], plots[[3]], plots[[4]],
  ncol =4, nrow = 1,
  common.legend = TRUE, legend = "right"
)

print(combined_plot)




