rm(list = ls())
library(CellChat)
library(Seurat)
library(tidyverse)
library(NMF)
library(ggalluvial)
library(patchwork)
library(ggplot2)
library(svglite)

scRNA <- readRDS("Mac_celltype.rds")
table(scRNA@meta.data$Sample)
Idents(scRNA) <- "celltype"
scRNA@commands$FindClusters
library(mindr)
(out <- capture.output(str(scRNA)))
out2 <- paste(out,collapse = "\n")
mm(gsub("\\.\\.@","#",gsub("\\.\\.","#",out2)),type ="text",root="Seurat")
DefaultAssay(scRNA) <- "RNA"
cellchat <- createCellChat(scRNA,group.by = "celltype")
cellchat
summary(cellchat)
str(cellchat)
(out <- capture.output(str(cellchat)))
out2 <- paste(out, collapse="\n")
mm(gsub("\\.\\.@","# ",gsub("\\.\\. ","#",out2)),type ="text")

levels(cellchat@idents)
groupsize <- as.numeric(table(cellchat@idents))
groupsize
cellchatDB <- CellChatDB.human
str(cellchatDB)
colnames(cellchatDB$interaction)
cellchatDB$interaction[1:4,1:4]
head(cellchatDB$geneInfo)
showDatabaseCategory(cellchatDB)
unique(cellchatDB$interaction$annotation)

cellchatDB.use <- subsetDB(cellchatDB,search = "Secreted Signaling")

#cellchatDB.use <- subsetDB(cellchatDB)
cellchat@DB <- cellchatDB.use

cellchat <- subsetData(cellchat)
future::plan("multiprocess",workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat,PPI.human)
cellchat <- computeCommunProb(cellchat,raw.use = F,population.size = T)
cellchat <- filterCommunication(cellchat,min.cells = 10)
df.net <- subsetCommunication(cellchat)
write.csv(df.net,"df.net.csv")
cellchat <- computeCommunProbPathway(cellchat)
df.netp <- subsetCommunication(cellchat,slot.name = "netP")
saveRDS(cellchat,"cellchat.rds")
cellchat@data@Dimnames[[1]] <- NULL
cellchat@data@Dimnames[[2]] <- NULL
cellchat@data@x <- 0
cellchat@data@i <- 0L
cellchat@data@p <- 0L
cellchat@data@Dim <- 0L

saveRDS(cellchat,"cellchat.rds")
