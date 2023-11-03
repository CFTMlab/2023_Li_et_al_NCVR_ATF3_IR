library(Seurat)
library(cowplot)
library(patchwork)
library(dplyr)
library(ggplot2)
chat1;chat2
Object.list <- list(chat1,chat2)
object.list <- lapply(X=Object.list, FUN=function(x){
  x<-NormalizeData(x)
  x<-FindVariableFeatures(x, selection.method = "vst", nfeature=2000)
})
integration <- function(object_list, dim_Num=20){
  dim_Num <- as.numeric(dim_Num)
  all_samples<-names(object_list)
  
  immune.anchors <- FindIntegrationAnchors(object.list = object_list, dims = 1:dim_Num)
  combined_Obj <- IntegrateData(anchorset = immune.anchors, dims = 1:dim_Num)
  DefaultAssay(combined_Obj) <- "integrated"
  
  combined_Obj <- ScaleData(combined_Obj, verbose = FALSE)
  combined_Obj <- RunPCA(combined_Obj, npcs = 20,verbose = FALSE)
  
  return (combined_Obj)
}

scRNA <- integration(Object.list)

gb1 <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(gb1), 10)
gbplotwt1 <- VariableFeaturePlot(gb1)
gbplotwt2 <- LabelPoints(plot = gbplotwt1, points = top10, repel = TRUE)
gbplotwt1
gbplotwt2
#ggsave("VariableFeatures.png", plot = gbplotwt1+gbplotwt2 , width = 15, height = 7)

all.genes <- rownames(scRNA)
scRNA <- ScaleData(scRNA, features = all.genes)
GetAssayData(scRNA,slot="counts",assay="RNA")          
GetAssayData(scRNA,slot="data",assay="RNA")
CaseMatch(c(cc.genes$s.genes,cc.genes$g2m.genes),VariableFeatures(scRNA))
g2m_genes = cc.genes$g2m.genes
g2m_genes = CaseMatch(search = g2m_genes, match = rownames(scRNA))
s_genes = cc.genes$s.genes
s_genes = CaseMatch(search = s_genes, match = rownames(scRNA))
scRNA<- CellCycleScoring(object=scRNA,  g2m.features=g2m_genes,  s.features=s_genes)
scRNA <- RunPCA(scRNA, features = c(s_genes, g2m_genes))
p <- DimPlot(scRNA, reduction = "pca", group.by = "Phase")
p

#PCA
scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA)) 
plot1 <- DimPlot(scRNA, reduction = "pca", group.by="orig.ident") 
plot2 <- ElbowPlot(scRNA, ndims=20, reduction="pca") 
plot1+plot2
#ggsave("PCA.png", plot=plot1 + plot2, width=8, height=6)
scRNA<- JackStraw(scRNA, num.replicate = 100)
scRNA <- ScoreJackStraw(scRNA, dims = 1:20)
JackStrawPlot(scRNA, dims = 1:20)
ElbowPlot(scRNA)

scRNA  <- FindNeighbors(scRNA , dims = 1:20)
scRNA  <- FindClusters(scRNA , resolution = 0.3)
scRNA <- RunUMAP(scRNA, dims = 1:20)
embed_umap <- Embeddings(scRNA, 'umap')  
#write.csv(embed_umap,'cluster1/embed_umap.csv') 
#group_by_cluster
DimPlot(scRNA, reduction = "umap",group.by = "integrated_snn_res.0.3", label=T) 
#ggsave("UMAP.png", plot = plot3, width = 8, height = 7)
#group_by_sample
DimPlot(scRNA, reduction = "umap", group.by='orig.ident')
#ggsave("UMAP_sample.png", plot = plot3 + plot4, width = 14, height = 7)
saveRDS(scRNA,"scRNA_celltype.rds")