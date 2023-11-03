library(Seurat)
library(cowplot)
library(patchwork)
library(dplyr)
library(ggplot2)
#cell number /// cell ratio
mac <- readRDS("Mac_celltype.rds")
DimPlot(mac,label = T)
Idents(mac) <- "celltype3"
as.data.frame(prop.table(table(Idents(mac), mac@meta.data[,"orig.ident"]), margin = 2))-> pdf -> td
write.csv(td,"cellratio_celltype3.csv")
library(tidyverse)
library(pheatmap)
data <- spread(td,key = "Var2",value = "Freq")
rownames(data) <- data[,1]
data <- data[,-1]
data <- round(data,2)
pheatmap(data,
         cluster_rows = F,
         cluster_cols = F,
         scale = "none",
         display_numbers = T,
         fontsize_number = 12,
         fontsize_row = 14,
         fontsize_col = 14,
         angle_col = 90,
         color = colorRampPalette(c("white","blue"))(50))
write.csv(td,"ratio.csv")
data <- data.frame(table(mac@meta.data$celltype,mac@meta.data$orig.ident))
write.csv(data,"number.csv")