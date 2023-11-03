library(Seurat)
library(cowplot)
library(patchwork)
library(dplyr)
library(ggplot2)
chat1.10x <- Read10X('IR0h_H/')
chat1 <- CreateSeuratObject(counts = chat1.10x, 
                            project = "IR0h_H", 
                            min.cells = 30, 
                            min.features = 500)
chat1[["percent.mito"]] <- PercentageFeatureSet(chat1, pattern = "^mt-")
VlnPlot(chat1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
chat1 <- subset(chat1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mito < 30)
plot1 <- FeatureScatter(chat1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1
chat1 <- NormalizeData(chat1, normalization.method = "LogNormalize", scale.factor = 2000)

chat2.10x <- Read10X('IR6h_H/')
chat2 <- CreateSeuratObject(counts = chat2.10x, 
                            project = "IR6h_H", 
                            min.cells = 30, 
                            min.features = 500)
chat1[["percent.mito"]] <- PercentageFeatureSet(chat1, pattern = "^mt-")
VlnPlot(chat1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
chat1 <- subset(chat1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mito < 30)
plot1 <- FeatureScatter(chat1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1
chat2 <- NormalizeData(chat2, normalization.method = "LogNormalize", scale.factor = 2000)