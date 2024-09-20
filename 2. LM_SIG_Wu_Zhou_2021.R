rm(list = ls()); gc()
options(stringsAsFactors = F)
setwd("G:/LUAD_scRNASeq/")

.libPaths(c("D:/R/R-4.3.2/Seurat_v4", .libPaths()))
unloadNamespace("infercnv")
library(Seurat)
library(GSVA)
library(clusterProfiler)
library(org.Hs.eg.db)

setwd("G:/LUAD_scRNASeq/merge_set_cancercell2022")
NSCLC_cc = readRDS("data.rds")
meta_cc = data.frame(NSCLC_cc@meta.data)

table(meta_cc$origin)
table(meta_cc$ann_fine)
table(meta_cc$assay)
table(meta_cc$dataset)

scRNA_list = NSCLC_cc[,NSCLC_cc@meta.data[["dataset"]]=="Wu_Zhou_2021"]
count = GetAssayData(scRNA_list,assay = "RNA",slot = "counts")
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
meta_data = scRNA_list@meta.data
all(rownames(meta_data) %in% colnames(seurat_obj))
meta_data <- meta_data[colnames(seurat_obj), ]
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, meta_data)

# read the list of LM-related genes
setwd("G:/LUAD_scRNASeq/")
Lactate_Metabolism_gs <- read.csv("lactate_gene.csv")
Lactate_Metabolism_gs = list(LM = Lactate_Metabolism_gs$sym)

rm(list = setdiff(ls(),c('seurat_obj',"Lactate_Metabolism_gs","meta_cc","NSCLC_cc","scRNA_list")))

sce <- seurat_obj
Idents(sce) <- "Celltype..malignancy."
sce <- sce[,sce@meta.data[["ann_fine"]]=="Cancer cells"]
counts <- sce@assays$RNA@data
dim(counts)
LM_gsva <- gsva(
  expr = counts, gset.idx.list = Lactate_Metabolism_gs, kcdf = "Gaussian",
  parallel.sz = 60
)
LM_gsva <- as.data.frame(t(LM_gsva))
LM_gsva_auc = as.data.frame(t(getAUC(LM_gsva)))
lmGenes <- data.frame(
  gene = rownames(counts),
  coef = NA, p = NA)
for (i in 1:nrow(counts)) {
  cor <- cor.test(counts[i, ], LM_gsva$LM, method = "spearman")
  lmGenes$coef[i] <- cor$estimate
  lmGenes$p[i] <- cor$p.value
}
lmGenes$p.adjust <- p.adjust(lmGenes$p, method = "BH")


sce <- scRNA_list
sce$group <- ifelse(
    scRNA_list@meta.data[["ann_fine"]] == "Cancer cells", "Malignant","control")
Idents(sce) <- "group"
future::plan("multisession", workers = 8)
DE <- FindMarkers(
      sce, ident.1 = "Malignant", group.by = "group", logfc.threshold = 0.25,
      min.pct = 0.1, base = exp(1)
      )
DE <- DE[DE$p_val_adj < 1e-05, ]
ids=bitr(rownames(DE),fromType="ENSEMBL",toType="SYMBOL",OrgDb=org.Hs.eg.db)
DE$symbol = ids$SYMBOL[match(rownames(DE),ids$ENSEMBL)]

LMx = lmGenes
LMy = DE
LMx <- LMx[LMx$coef > 0 & LMx$p.adjust < 1e-05, ]
LMy <- LMy$symbol[LMy$avg_logFC >= 0.25]
LMy <- LMy[!grepl("^RP[SL]", LMy, ignore.case = F)]
LMn <- LMx[LMx$gene %in% LMy, ]

save(LM_gsva,lmGenes,file = "cc_Wu_Zhou_2021_LMx_LMy.RData")
save(LMn,file = "cc_Wu_Zhou_2021_LMn.RData")





