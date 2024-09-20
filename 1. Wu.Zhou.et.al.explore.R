rm(list = ls()); gc()
options(stringsAsFactors = F)
setwd("G:/LUAD_scRNASeq/")

.libPaths(c("D:/R/R-4.3.2/Seurat_v4", .libPaths()))
unloadNamespace("infercnv")
library(Seurat)
library(GSVA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

setwd("G:/LUAD_scRNASeq/merge_set_cancercell2022")
NSCLC_cc = readRDS("data.rds")
meta_cc = data.frame(NSCLC_cc@meta.data)
table(meta_cc$origin)
table(meta_cc$ann_fine)
table(meta_cc$assay)
table(meta_cc$dataset)
scRNA = NSCLC_cc[,NSCLC_cc@meta.data[["dataset"]]=="Wu_Zhou_2021"]
rm(NSCLC_cc)
rm(meta_cc)
meta = scRNA@meta.data
count = GetAssayData(scRNA,assay = "RNA",slot = "counts")
ensembl=rownames(count)
ids=bitr(ensembl,fromType="ENSEMBL",toType="SYMBOL",OrgDb=org.Hs.eg.db)
rownames(count) = ids$SYMBOL[match(rownames(count),ids$ENSEMBL)]
count=count[!is.na(rownames(count)),]
du = unique(rownames(count)[duplicated(rownames(count))])
aggregated_expr <- matrix(0, nrow = length(du), ncol = ncol(count))
rownames(aggregated_expr) <- du
colnames(aggregated_expr) <- colnames(count)
for (gene in du) {
  aggregated_expr[gene, ] <- colSums(count[rownames(count) == gene, ])
}
non_duplicated_expr <- count[!rownames(count) %in% du, ]
final_expr_matrix <- rbind(non_duplicated_expr, aggregated_expr)
count = final_expr_matrix


seurat_obj <- CreateSeuratObject(counts = count, project = "MyProject")
meta_data = scRNA@meta.data
all(rownames(meta_data) %in% colnames(seurat_obj))
meta_data <- meta_data[colnames(seurat_obj), ]
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, meta_data)
save(seurat_obj,file = "Wu.Zhou.et.al.RData")  

rm(list=ls())
gc()

setwd("G:/LUAD_scRNASeq/merge_set_cancercell2022")
load("Wu.Zhou.et.al.RData")
load("G:/LUAD_scRNASeq/cc_Wu_Zhou_2021_LMn.RData")
LMn = intersect(LMn$gene,rownames(seurat_obj))
LMn = list(LM = LMn)

meta = seurat_obj@meta.data
cell_meta = data.frame(table(meta[,"ann_coarse"]))
colnames(cell_meta) = c("cell type","cell number")
cell_meta$ratio = cell_meta$`cell number`/sum(cell_meta$`cell number`)
write.csv(cell_meta,"G:/LUAD_scRNASeq/Fig_and_chart/Fig1A.csv",row.names = F)

library(irGSEA)
library(AUCell)
library(UCell)
library(GSEABase)
genes <- LMn[[1]]
geneSets <- GeneSet(genes, setName="LM")
seurat_obj <- irGSEA.score(object = seurat_obj,assay = "RNA", 
                             slot = "data", seeds = 123, ncores = 1,
                             min.cells = 3, min.feature = 0,
                             custom = T, geneset = LMn, msigdb = F, 
                             geneid = "symbol",
                             method = "UCell",
                             aucell.MaxRank = NULL, ucell.MaxRank = NULL, 
                             kcdf = 'Gaussian')
sce.all = seurat_obj
sce.all[["percent.mt"]] <- PercentageFeatureSet(sce.all, pattern = "^MT-")
sce.all[["percent.rp"]] <- PercentageFeatureSet(sce.all, pattern = "^RP[SL]")
VlnPlot(sce.all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 2)
plot1 <- FeatureScatter(sce.all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sce.all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(sce.all, feature1 = "nCount_RNA", feature2 = "percent.rp")
plot1 + plot2 + plot3

sce.all[["HKG.score"]] = colSums(sce.all[c("ACTB","GAPDH","MALAT1"),])
sce.all.qc <- subset(sce.all, subset = nFeature_RNA > 500 & percent.mt < 20 & percent.rp < 50 & HKG.score >= 1)

sce <- SCTransform(sce.all)
sce <- RunPCA(sce, npcs=50, verbose=FALSE)
sce <- RunUMAP(sce, dims = 1:50, verbose = FALSE)
DimPlot(sce,reduction = "pca")
DimPlot(sce,reduction = "umap",group.by = "cell_type_predicted") + 
  ggtitle("Wu Zhou et al cell type")

gc()

sce@meta.data[["sample"]] = droplevels(sce@meta.data[["sample"]])
sce@meta.data[["donor_id"]] = droplevels(sce@meta.data[["donor_id"]])
sce@meta.data[["ann_fine"]] = droplevels(sce@meta.data[["ann_fine"]])
sce@meta.data[["cell_type_predicted"]] = droplevels(sce@meta.data[["cell_type_predicted"]])
sce@meta.data[["ann_coarse"]] = droplevels(sce@meta.data[["ann_coarse"]])
sce@meta.data[["cell_type_tumor"]] = droplevels(sce@meta.data[["cell_type_tumor"]])
sce@meta.data[["cell_type_tumor"]] = droplevels(sce@meta.data[["cell_type_tumor"]])

DimPlot(sce,reduction = "umap",group.by = "cell_type_predicted") + 
  ggtitle("Wu Zhou et al cell type")


library(harmony)
sce <- RunHarmony(sce, group.by.vars="sample", max.iter.harmony =
                    20,assay.use = "SCT")
ElbowPlot(sce, ndims=50, reduction="pca")

pc.num=1:50
sce <- RunTSNE(sce, reduction="harmony", dims=pc.num) %>%
  RunUMAP(reduction="harmony", dims=pc.num)
sce <- FindNeighbors(sce, dims = pc.num)
sce <- FindClusters(sce, resolution = 0.6)
table(sce@meta.data$seurat_clusters)

DimPlot(sce, reduction = "umap",label = T,group.by = "ann_coarse",label.size = 4) + 
  NoLegend() + ggtitle("Wu Zhou et al cell type")+
  theme(plot.title = element_text(size = 18))
DimPlot(sce, reduction = "umap",label = T,group.by = "ann_fine") + NoLegend()


DimPlot(sce, reduction = "umap",label = T,group.by = "ann_coarse") + NoLegend() + ggtitle("Wu Zhou et al cell type")
UCell = as.data.frame(sce@assays[["UCell"]]@counts)
UCell = as.data.frame(t(UCell[,Cells(seurat_obj)]))
colnames(UCell) = "LM_UCell"
sce@meta.data$LM_UCell = UCell$LM_UCell
FeaturePlot(sce,features = "LM_UCell") + ggtitle("lactate metabolism strength")

UCell_meta = sce@meta.data
UCell_meta = UCell_meta[,c("ann_coarse","LM_UCell")]
median_UCell = UCell_meta %>% group_by(ann_coarse) %>% summarise(median_score = median(LM_UCell))
pair_p = pairwise.t.test(UCell_meta$LM_UCell,UCell_meta$ann_coarse,p.adjust.method = "none")
pair_p = pair_p[["p.value"]]

write.csv(median_UCell,"G:/LUAD_scRNASeq/Fig_and_chart/Fig1C_median.csv",row.names = F)
write.csv(pair_p,"G:/LUAD_scRNASeq/Fig_and_chart/Fig1C_pairwise_p.csv",row.names = T)
meta_sce = sce@meta.data

VlnPlot(sce,features = "LM_UCell",group.by = "ann_fine",pt.size = 0)

library(dplyr)
mean_scores <- meta_sce %>%
  group_by(ann_coarse) %>%   
  summarise(mean_score = mean(LM_UCell, na.rm = TRUE))

sce = SetIdent(sce,value = "ann_coarse")
markers = FindAllMarkers(sce,group.by = "ann_coarse",only.pos = T)
LMn = LMn[[1]]
markers_LM = markers[markers$gene %in% LMn & markers$cluster != "Epithelial cell",]
feature_gene = c("IFI27","IGFBP3","TM4SF1","ENAH")
FeaturePlot(sce,features = feature_gene)

write.csv(markers_LM,"LM_gene_in_non_magliant_cells.csv",row.names = T)


sce_no_epi = sce[,sce@meta.data$ann_coarse != "Epithelial cell"]
sce_no_epi <- SCTransform(sce_no_epi)
sce_no_epi <- RunPCA(sce_no_epi, npcs=50, verbose=FALSE)
sce_no_epi <- RunUMAP(sce_no_epi, dims = 1:50, verbose = FALSE)
DimPlot(sce_no_epi,reduction = "umap",group.by = "ann_coarse",label = T,label.size = 4.5) + NoLegend() +
  ggtitle("Wu Zhou et al cell type without magliant cells") + 
  theme(plot.title = element_text(size = 18))

feature_gene = c("IFI27","IGFBP3","TM4SF1","ENAH")
FeaturePlot(sce_no_epi,features = feature_gene)
VlnPlot(sce_no_epi,features = feature_gene,group.by = "ann_coarse",pt.size = 0.1)

markers_LM = markers_LM[feature_gene,]
sce_no_epi_count = as.data.frame(t(sce_no_epi@assays$RNA$counts))
sce_no_epi_count = sce_no_epi_count[,feature_gene]
sce_no_epi_count$cell_type = sce_no_epi@meta.data$ann_coarse

res = matrix(0,nrow = length(feature_gene),ncol = length(unique(sce_no_epi_count$cell_type)),
             dimnames = list(feature_gene,unique(sce_no_epi_count$cell_type)))
for(gene in feature_gene){
  for(cell_type in unique(sce_no_epi_count$cell_type)){
    df = sce_no_epi_count[sce_no_epi_count$cell_type == cell_type,gene]
    por = sum(df>0)/length(df)
    res[gene,cell_type] = por
  }
}

all(markers_LM$gene == rownames(res))
rownames(markers_LM) = markers_LM$gene
colnames(res) = paste0("pct.in " ,colnames(res))
markers_LM = cbind(markers_LM,res)
write.csv(pair_p,"G:/LUAD_scRNASeq/Fig_and_chart/Fig1D_4_gene_marker.csv",row.names = T)
