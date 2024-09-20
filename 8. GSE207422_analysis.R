.libPaths(c("D:/R/R-4.3.2/Seurat_v4", .libPaths()))
unloadNamespace("infercnv")
library(Seurat)

setwd("G:/LUAD_scRNASeq/GSE207422/")
#整合cell-type注释
rm(list=ls())
gc()
load("G:/MR_mGWAS/GSE207422/GSE207422_Mye_subcluster.RData")
load("G:/MR_mGWAS/GSE207422/GSE207422_TNK_subcluster.RData")
load("G:/MR_mGWAS/GSE207422/GSE207422_qc.RData")
cell_meta = sce@meta.data
M_cell_meta = M_cell@meta.data
T_cell_meta = T_cell@meta.data
TM_cell_meta = rbind(T_cell_meta[,c("orig.ident","nCount_RNA","nFeature_RNA","type")],M_cell_meta[,c("orig.ident","nCount_RNA","nFeature_RNA","type")])
cell_meta[rownames(TM_cell_meta),"type"] = TM_cell_meta$type
cell_meta$type[is.na(cell_meta$type)] = cell_meta$singleR_hpca[is.na(cell_meta$type)]
cell_meta$type[cell_meta$customclassif == "Naive B cells"] = "Naive B cells"
cell_meta$type[cell_meta$customclassif == "Plasma B cells"] = "Plasma B cells"
all(colnames(sce) == rownames(cell_meta))
sce$final_celltype = cell_meta$type
DimPlot(sce, reduction='umap', group.by="final_celltype",label = T,label.size = 4,repel = TRUE)


#免疫基因cell-type分布:"CCT6A","KRT8","LYPD3","HMGA1","PSMB7
select_genes = c("CCT6A","KRT8","LYPD3","HMGA1","PSMB7")
DotPlot(sce, features = select_genes,group.by  = "singleR_hpca")+ 
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) #颜色

DotPlot(sce, features = select_genes,group.by  = "final_celltype")+ 
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) #颜色
FeaturePlot(sce, features = select_genes)
VlnPlot(sce, features = select_genes, ncol = 2,group.by = "final_celltype")

#输出下raw分群和整体分群的细胞数量与比例
Fig4A = cell_meta[,c("singleR_hpca","type")]
Fig4A1 = data.frame(table(Fig4A$singleR_hpca))
Fig4A2 = data.frame(table(Fig4A$type))
Fig4A1$ratio = as.numeric(Fig4A1$Freq)/sum(as.numeric(Fig4A1$Freq))
Fig4A2$ratio = as.numeric(Fig4A2$Freq)/sum(as.numeric(Fig4A2$Freq))
colnames(Fig4A1) = c("cell_type","cell_num","cell_ratio")
colnames(Fig4A2) = c("cell_type","cell_num","cell_ratio")

#存到2个excel表中
library(openxlsx)
wb = createWorkbook()
addWorksheet(wb,"celltype_raw")
writeData(wb,"celltype_raw",Fig4A1)
addWorksheet(wb,"celltype_fined")
writeData(wb,"celltype_fined",Fig4A2)
saveWorkbook(wb,"G:/LUAD_scRNASeq/Fig_and_chart/Fig4A.xlsx",overwrite = T)

#计算5-gene的非零比例
marker_mat = sce@assays[["RNA"]]@counts
marker_mat = as.matrix(marker_mat[c("CCT6A","KRT8","LYPD3","HMGA1","PSMB7"),])
marker_df = as.data.frame(t(marker_mat))
marker_df$cell_type = cell_meta$type

res = marker_df %>% group_by(cell_type) %>% 
  summarise(across(everything(),~sum(.>0)/n(),.names="prop_{col}"))

write.csv(res,"G:/LUAD_scRNASeq/Fig_and_chart/Fig4B.csv",row.names = F)

#计算平均表达量(log2CPM)与PR、SD和浸润比例相关性
pseudo_sce = AggregateExpression(sce,assays = "SCT",group.by = "orig.ident",return.seurat = T)
pseudo_sce = as.data.frame(pseudo_sce@assays[["SCT"]]@counts)
pseudo_sce = apply(pseudo_sce ,2, function(x) { x/sum(x)*1e6 })
cell_ratio = data.frame(table(cell_meta[,c("orig.ident",'type')]))
cell_sum = data.frame(table(cell_meta$orig.ident))
cell_ratio$sum = cell_sum$Freq[match(cell_ratio$orig.ident,cell_sum$Var1)]
cell_ratio$ratio = cell_ratio$Freq/cell_ratio$sum
cell_ratio$CCT6A = pseudo_sce["CCT6A",cell_ratio$orig.ident]
cell_ratio$KRT8 = pseudo_sce["KRT8",cell_ratio$orig.ident]
cell_ratio$LYPD3 = pseudo_sce["LYPD3",cell_ratio$orig.ident]
cell_ratio$HMGA1 = pseudo_sce["HMGA1",cell_ratio$orig.ident]
cell_ratio$PSMB7 = pseudo_sce["PSMB7",cell_ratio$orig.ident]
cell_ratio$RECIST = cell_meta$RECIST[match(cell_ratio$orig.ident,cell_meta$orig.ident)]

library(tidyr)
library(dplyr)
library(psych)
library(ggplot2)
library(ggsignif)
library(reshape2)
library(pheatmap)
library(ggpubr)

cell_type = names(table(cell_ratio$type))
cor_matrix <- matrix(nrow = length(select_genes), ncol = length(cell_type))
p_matrix <- matrix(nrow = length(select_genes), ncol = length(cell_type))

# 计算相关系数和 p-value
for (i in 1:length(select_genes)) {
  for (j in 1:length(cell_type)) {
    df = cell_ratio[cell_ratio$type == cell_type[j],]
    test_result <- cor.test(df[, colnames(df)==select_genes[i]], df[, "ratio"])
    cor_matrix[i, j] <- test_result$estimate
    p_matrix[i, j] <- test_result$p.value
  }
}

colnames(cor_matrix) = cell_type
rownames(cor_matrix) = select_genes
colnames(p_matrix) = cell_type
rownames(p_matrix) = select_genes


#画图时不运行这些
cor_matrix = data.frame(cor_matrix)
p_matrix = data.frame(p_matrix)
cor_matrix = cbind(rownames(cor_matrix),cor_matrix)
p_matrix = cbind(rownames(p_matrix),p_matrix)
colnames(cor_matrix)[1] = "LM_gene"
colnames(p_matrix)[1] = "LM_gene"

#存到2个excel表中
library(openxlsx)
wb = createWorkbook()
addWorksheet(wb,"pcc")
writeData(wb,"pcc",data.frame(cor_matrix))
addWorksheet(wb,"pcc_pvalue")
writeData(wb,"pcc_pvalue",data.frame(p_matrix))
saveWorkbook(wb,"G:/LUAD_scRNASeq/Fig_and_chart/Fig4C.xlsx",overwrite = T)

# 转置矩阵以匹配要求的维度
cor_matrix = cor_matrix[,-1] #把添加的第一列基因名去除
p_matrix = p_matrix[,-1]
cor_matrix[is.na(cor_matrix)] = 0
p_matrix[is.na(p_matrix)] = 1
# 使用 pheatmap 绘制热图
significance_marks<- matrix(0, nrow = 5, ncol = 30)# 初始化空字符串矩阵
significance_marks <- ifelse(p_matrix < 0.05, "*", "")

significance_marks[p_matrix < 0.001] <- "***"
significance_marks[p_matrix >= 0.001 & p_matrix < 0.01] <- "**"
library(ComplexHeatmap)
pheatmap(cor_matrix,display_numbers = significance_marks,cluster_rows = FALSE,
         cluster_cols = FALSE,show_rownames = TRUE,show_colnames = TRUE,
         color = colorRampPalette(c("#2C7BB6", "#FFFFBF", "#D7191C"))(100),
         main = "LM risk gene PCC with TME cells"
)

# 设置输出的PDF文件名
pdf("Fig4B-5-LM-risk-gene-PCC_TME_cells.pdf", width = 7, height = 5)
# 绘制多个图
pheatmap(cor_matrix,display_numbers = significance_marks,cluster_rows = FALSE,
         cluster_cols = FALSE,show_rownames = TRUE,show_colnames = TRUE,
         color = colorRampPalette(c("#2C7BB6", "#FFFFBF", "#D7191C"))(100),
         main = "LM risk gene PCC with TME cells"
)
# 关闭图形设备，完成输出
dev.off()

#PD/SR和浸润、表达的boxplot
ggplot(cell_ratio,aes(x = type,y=ratio,fill=RECIST)) + 
  geom_boxplot() + theme_bw() + theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  labs(title = "immune porp vs aPD1 response",x="cellType",y="porp",fill = "RECIST") +
  stat_compare_means(aes(group = RECIST),label = "p.signif")

df = pivot_longer(cell_ratio,cols=select_genes,names_to = "genes",values_to = "exp")


# 设置输出的PDF文件名
pdf("Fig4C-5-LM-risk-gene-vs-aPD1-response.pdf", width = 7, height = 5)
# 绘制多个图
ggplot(df,aes(x = genes,y=exp,fill=RECIST)) + 
  geom_boxplot(outlier.shape = NA) + theme_bw() + theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  scale_y_continuous(limits =c(0,300)) + 
  labs(title = "LM risk gene vs aPD1 response",x="risk_gene",y="exp",fill = "RECIST")+
  stat_compare_means(aes(group = RECIST),label = "p.signif")
dev.off()

#读取bulk RNA-seq，转化counts，计算risk score
library(readxl)
library(purrr)
bulk = read.table("G:/MR_mGWAS/GSE207422/GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt",sep="\t",header=T)
bulk_meta = read_excel("G:/MR_mGWAS/GSE207422/GSE207422_NSCLC_bulk_RNAseq_metadata.xlsx")
rownames(bulk) = bulk$Gene
bulk = bulk[,-1]
bulk_meta = na.omit(bulk_meta)
all(colnames(bulk) == bulk_meta$Sample)
bulk_meta$CCT6A = t(bulk["CCT6A",])
bulk_meta$KRT8 = t(bulk["KRT8",])
bulk_meta$LYPD3 = t(bulk["LYPD3",])
bulk_meta$HMGA1 =t(bulk["HMGA1",])
bulk_meta$PSMB7 =t(bulk["PSMB7",])

colnames(bulk_meta)[13:17] = c("CCT6A","KRT8","LYPD3","HMGA1","PSMB7")

df = pivot_longer(bulk_meta,cols=select_genes,names_to = "genes",values_to = "exp")
df = df[df$RECIST %in% c("SD","PR"),]
t=ggplot(df,aes(x = genes,y=exp,fill=RECIST)) + 
  geom_boxplot(outlier.shape = NA) + theme_bw() + theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  scale_y_continuous(limits =c(0,10)) + 
  labs(title = "risk gene vs aPD1 response in Bulk RNA-Seq",x="risk_gene",y="exp",fill = "RECIST")+
  stat_compare_means(aes(group = RECIST),label = "p.signif")

# 设置输出的PDF文件名
pdf("Fig4C-5-LM-risk-gene-vs-aPD1-response-bulk.pdf", width = 7, height = 5)
# 绘制多个图
t
# 关闭图形设备，完成输出
dev.off()


#存储图4D的两个boxplot的两个表
Fig4DA = cell_ratio[,6:11]
median_4DA = Fig4DA %>% group_by(RECIST) %>% summarise(across(1:5,median))
p_4DA = map_dfr(1:5,function(i){
  p.value = wilcox.test(Fig4DA[[i]]~Fig4DA$RECIST)$p.value
  data.frame(Score = colnames(Fig4DA)[i],p_value = p.value)
})
Fig4DB = bulk_meta[bulk_meta$RECIST %in% c("SD","PR"),12:17]
median_4DB = Fig4DB %>% group_by(RECIST) %>% summarise(across(1:5,median))
p_4DB = map_dfr(2:6,function(i){
  p.value = wilcox.test(Fig4DB[[i]]~Fig4DB$RECIST)$p.value
  data.frame(Score = colnames(Fig4DB)[i],p_value = p.value)
})

blank = data.frame(matrix(0,nrow = 1,ncol = 6))
colnames(blank) = colnames(median_4DA)
median_4DA = rbind(median_4DA,blank)
median_4DB = rbind(median_4DB,blank)

median_4DA[3,2:6] = t(p_4DA$p_value)
median_4DB[3,2:6] = t(p_4DB$p_value)
median_4DA[3,1] = "pvalue"
median_4DB[3,1] = "pvalue"


#存到2个excel表中
library(openxlsx)
wb = createWorkbook()
addWorksheet(wb,"risk gene in scRNA-Seq")
writeData(wb,"risk gene in scRNA-Seq",median_4DA)
addWorksheet(wb,"risk gene in bulk RNA-Seq")
writeData(wb,"risk gene in bulk RNA-Seq",median_4DB)
saveWorkbook(wb,"G:/LUAD_scRNASeq/Fig_and_chart/Fig4D.xlsx",overwrite = T)


#为防止重复，做个新的GSE207422的散点图
load("G:/MR_mGWAS/GSE207422/GSE207422_qc.RData")
sce <- SCTransform(sce)
Idents(sce) = sce$orig.ident
sce <- RunPCA(sce, npcs=50, verbose=FALSE)
sce <- RunUMAP(sce, dims = 1:50, verbose = FALSE)
DimPlot(sce,reduction = "pca")
DimPlot(sce,reduction = "umap")
rm("sce.all.qc")
gc()

#这里多个单细胞样本整合的方法为harmony，也可以用ScVI
library(harmony)
sce <- RunHarmony(sce, group.by.vars="orig.ident", max.iter.harmony =
                    20,assay.use = "SCT")
#选择PC数量的方法：Elbowplot
ElbowPlot(sce, ndims=50, reduction="pca")

pc.num=1:40
sce <- RunTSNE(sce, reduction="harmony", dims=pc.num) %>%
  RunUMAP(reduction="harmony", dims=pc.num)
sce <- FindNeighbors(sce, dims = pc.num)
#res不同，做个不一样的散点图，还可以读取res=0.6时的细胞亚群
sce <- FindClusters(sce, resolution = 0.6)
table(sce@meta.data$seurat_clusters)
DimPlot(sce, reduction = "tsne",label = F,group.by = "orig.ident")
DimPlot(sce, reduction = "umap",label = F,group.by = "orig.ident")


load("G:/MR_mGWAS/GSE207422/GSE207422_Mye_subcluster.RData")
load("G:/MR_mGWAS/GSE207422/GSE207422_TNK_subcluster.RData")
cell_meta = sce@meta.data
M_cell_meta = M_cell@meta.data
T_cell_meta = T_cell@meta.data
TM_cell_meta = rbind(T_cell_meta[,c("orig.ident","nCount_RNA","nFeature_RNA","type")],M_cell_meta[,c("orig.ident","nCount_RNA","nFeature_RNA","type")])
cell_meta[rownames(TM_cell_meta),"type"] = TM_cell_meta$type
cell_meta$type[is.na(cell_meta$type)] = cell_meta$singleR_hpca[is.na(cell_meta$type)]
cell_meta$type[cell_meta$customclassif == "Naive B cells"] = "Naive B cells"
cell_meta$type[cell_meta$customclassif == "Plasma B cells"] = "Plasma B cells"
all(colnames(sce) == rownames(cell_meta))
sce$final_celltype = cell_meta$type
library(ggplot2)
DimPlot(sce, reduction='umap', group.by="singleR_hpca",label = T,label.size = 3.3,repel = TRUE) + NoLegend() + ggtitle("GSE207422 raw cell type")
DimPlot(sce, reduction='umap', group.by="final_celltype",label = T,label.size = 3.3,repel = TRUE) + NoLegend() + ggtitle("GSE207422 refined cell type")


pdf("Fig4A-GSE20-dimplot.pdf", width = 9, height = 7)
# 绘制多个图
DimPlot(sce, reduction='umap', group.by="final_celltype",label = T,label.size = 3.3,repel = TRUE) + NoLegend() + ggtitle("GSE207422 refined cell type")
dev.off()


meta = sce@meta.data
meta = meta[,c("singleR_hpca","seurat_clusters")]
#整体、T细胞、单核细胞的marker dotplot
library(ggplot2)
select_genes <- c("TRAC","CD3G","CD3E","CD3D",   #T cell
                  "KLRD1","GNLY","NKG7","NCAM1", #NK
                  "LYZ","FCGR3A","MARCO","CD68",  #Macrophage
                  "CLEC4E","CSF3R","CXCL8","FCGR3B", #Neutrophils
                  "IGHA2","IGHG3","IGHM","CD79A",  #B cell
                  "KRT18","CDH1","KRT19","EPCAM", #Epithelial
                  "THY1","COL1A2","COL1A1","DCN", #Fibroblast
                  "RAMP2","FLT1","CLDN5","PECAM1", #Endothelial
                  "GATA2","MS4A2","KIT",          #肥大细胞
                  "CD83","CD141","CD1A","CD1C",     #DC
                  "OLIG2","CNP","FA2H","MBP"    #少突胶质细胞
)

#顺序：T NK Mono/Macro Neutrophils B Epi Fibo Endo Mast DC
DotPlot(sce, features = select_genes )+ coord_flip() +
  scale_y_discrete(limits = c(
    '0','1','2','7','13','17','20',
    '9',
    '3','10','11','12','18','27',
    '4','8','23',
    '5','15','21','22','28','29','33',
    '6', '14', '16', '19', '26','30',
    '25',
    '31',
    '24',
    '32'
  ))


