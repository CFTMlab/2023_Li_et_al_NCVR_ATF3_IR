library(Seurat)
library(cowplot)
library(patchwork)
library(dplyr)
library(ggplot2)
mac <- readRDS("Mac_celltype.rds")
mac@meta.data$celltype3 <- factor(mac@meta.data$seurat_clusters,
                                  labels = c("BLT1_Mac","Trem2_Mac","MHCII_Mac","Ifit2_mac","Lyve1_Mac","Ly6c2_Mon","MHCII_Mac","Ki67_Mac","S100a9_Mon","Ly6c2_Mon","Ki67_Mac","Ki67_Mac"))
DimPlot(mac,label = T,split.by =  "orig.ident")
head(mac@meta.data)
gene <- c("Mertk","Ccr2")
VlnPlot(mac,features = c("rna_Mertk","rna_Ccr2"),pt.size = 0.1)
select_top_n<-function(scores,n_top){
  d <- data.frame(
    x   = data.table::copy(scores),
    indice=seq(1,length(scores)))
  
  data.table::setDT(d)
  data.table::setorder(d,-x)
  n_top_indice<-d$indice[1:n_top]
  return(n_top_indice)
}
cosg<-function(
    object,
    groups='all',
    assay='RNA',
    slot='data',
    mu=1,
    n_genes_user=100
){
  
  ### Obtain the cellxgene data
  genexcell <- Seurat::GetAssayData(object = object[[assay]], slot = slot)
  
  if (groups == 'all'){
    group_info <- Seurat::Idents(object = object)
  }else{ 
    object <- subset(x = object, idents = groups)
    group_info <- Seurat::Idents(object = object)
  }
  
  
  ### unique groups
  groups_order=sort(unique(group_info))
  n_cluster=length(groups_order)
  
  if (n_cluster == 1){
    stop('Cannot perform marker gene identification on a single cluster.')}
  
  
  n_cell=ncol(genexcell)
  n_gene=nrow(genexcell)
  gene_name=rownames(genexcell)
  
  ### If sepcifying too many genes to return
  if (n_genes_user>n_gene){
    n_genes_user=n_gene
  }
  
  
  cluster_mat=matrix(0,nrow =n_cluster,ncol = n_cell)
  
  order_i=1
  ### Set gene lambda and gene omega
  for (group_i in groups_order){
    idx_i=group_info==group_i 
    cluster_mat[order_i,idx_i]=1
    order_i=order_i+1
  }
  
  
  cluster_mat_sparse=as(cluster_mat, "dgCMatrix")
  ### Calculate the cosine similarity
  cosine_sim=proxyC::simil(genexcell,cluster_mat_sparse, method = "cosine",drop0=TRUE)
  
  pos_nonzero = cosine_sim != 0
  pos_nonzero=which(as.matrix(pos_nonzero),arr.ind = TRUE)
  
  #### Second-stage
  genexlambda=cosine_sim*cosine_sim
  e_power2_sum=Matrix::rowSums(genexlambda)
  
  if (mu==1){
    genexlambda[pos_nonzero]=genexlambda[pos_nonzero]/(replicate(ncol(genexlambda),e_power2_sum)[as.matrix(pos_nonzero)])
  }else{
    genexlambda[pos_nonzero]=genexlambda[pos_nonzero]/((
      (1-mu)*genexlambda[pos_nonzero] + mu * (replicate(ncol(genexlambda),e_power2_sum))
    )[as.matrix(pos_nonzero)])
  }
  
  genexlambda=genexlambda*cosine_sim
  
  rank_stats_names=data.frame(matrix(matrix(), n_genes_user, length(groups_order),
                                     dimnames=list(seq(1,n_genes_user), groups_order)),
                              stringsAsFactors=F)
  rank_stats_scores=data.frame(matrix(matrix(), n_genes_user, length(groups_order),
                                      dimnames=list(seq(1,n_genes_user), groups_order)),
                               stringsAsFactors=F)
  
  order_i=1
  ### Set gene lambda and gene omega
  for (group_i in groups_order){
    idx_i=group_info==group_i 
    scores=genexlambda[,order_i]
    global_indices = select_top_n(scores, n_genes_user)
    rank_stats_names[,order_i]=gene_name[global_indices]
    rank_stats_scores[,order_i]=scores[global_indices]
    
    ### save the group names
    order_i=order_i+1
  }
  
  colnames(rank_stats_names) <- groups_order
  colnames(rank_stats_scores) <- groups_order
  
  ###
  ranks_stats=list(
    names=rank_stats_names,
    scores=rank_stats_scores
    
  )
  ### return
  return(ranks_stats)
}
Idents(mac) <- "celltype3"
marker_cosg <- cosg(
  mac,
  groups='all',
  assay='RNA',
  mu=1,
  n_genes_user=100)
write.csv(marker_cosg$names,"marker_celltype3_cosg.csv")
gene <- c(marker_cosg$names$Trem2_Mac[1:10],
          marker_cosg$names$MHCII_Mac[1:10],
          marker_cosg$names$Lyve1_Mac[1:10],
          marker_cosg$names$BLT1_Mac[1:10],
          marker_cosg$names$Ly6c2_mono[1:10],
          marker_cosg$names$S100a9_mono[1:10],
          marker_cosg$names$Ki67_Mac[1:10])
DefaultAssay(mac) <- "integrated"
DoHeatmap(mac,features = gene,group.by = "celltype3") + scale_fill_gradientn(colors = c("#0099CC","white","#880000"))

rawdata <- data.frame(mac@assays$RNA@counts)
expre <- rawdata[gene,]
expre <- data.frame(t(expre))
head(expre)
expre$celltype1 <- ifelse(expre$Mertk == 0,"Mertk-","Mertk+")
expre$celltype2 <- ifelse(expre$Ccr2 == 0,"Ccr2-","Ccr2+")
expre$celltype <- paste0(expre$celltype1,expre$celltype2)
mac@meta.data$celltype_new2 <- expre$celltype