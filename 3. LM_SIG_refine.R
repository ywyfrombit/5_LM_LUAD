#从366个初次选出的基因中进一步去除：分泌蛋白、在其他微环境细胞中也是标志物的基因
library(readxl)
HPCA_serect <- read_excel("C:/Users/ywy40/Downloads/HPCA_serect_pro.xlsx")
View(HPCA_serect)
class = HPCA_serect$`Protein class`
class = strsplit(class,split = ",")
View(class)
class = unique(unlist(class))
class
HPCA_serect = HPCA_serect[grepl("Predicted secreted proteins",HPCA_serect$`Protein class`),]
load("G:/LUAD_scRNASeq/cc_Wu_Zhou_2021_LMn.RData")
intersect(LMn$gene,HPCA_serect$Gene)
LMn = LMn[!LMn$gene %in% HPCA_serect$Gene,]
View(LMn)


#其他单细胞的标志物
markers = read.csv("G:/LUAD_scRNASeq/merge_set_cancercell2022/LM_gene_in_non_magliant_cells.csv",header = T)
markers = markers$gene[markers$cluster %in% c("B cell","Neutrophils","NK cell","pDC","T cell")]
LMn = LMn[!LMn$gene %in% markers,]
write.csv(LMn,"G:/LUAD_scRNASeq/Fig_and_chart/LM_gene_LUAD.csv",row.names = F)
