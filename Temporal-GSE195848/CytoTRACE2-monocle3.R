.libPaths()
.libPaths(c('~/R/x86_64-pc-linux-gnu-library/4.4/'))
getwd()
setwd("/home/data/t060202/Spatial-GSVA-duijian/immune landscape during GBM progression/")

library(tidyverse)
library(COSG)
library(harmony)
library(ggsci)
library(dplyr) 
library(future)
library(Seurat)
library(clustree)
library(cowplot)
library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)
library(CytoTRACE2) 
library(Seurat)
library(paletteer)
library(BiocParallel)
register(MulticoreParam(workers = 8, progressbar = TRUE))
options(future.globals.maxSize = 100* 1024^3)

sub_data <- readRDS('meyloid_SeuratObj.RDS')

# check一下
DimPlot(sub_data,pt.size = 0.8,group.by = "cell_type",label = T)

dir.create("./13-monocle")
setwd("./13-monocle")
#2.数据预处理
cytotrace2_res <- cytotrace2(sub_data, #seurat对象
                             is_seurat = TRUE, 
                             slot_type = "counts", #counts和data都可以
                             species = 'mouse')#物种要选择human，默认是

#3.可视化
annotation <- data.frame(phenotype = sub_data@meta.data$cell_type) %>% 
  set_rownames(., colnames(sub_data))

# plotting-一次性生成多个图，然后储存在一个list，用$查看即可
plots <- plotData(cytotrace2_result = cytotrace2_res, 
                  annotation = annotation, 
                  is_seurat = TRUE)


plots$CytoTRACE2_UMAP
ggsave("CytoTRACE2_UMAP.pdf",  width = 9, height = 7, dpi = 300)
plots$CytoTRACE2_Potency_UMAP
ggsave("CytoTRACE2_Potency_UMAP.pdf",  width = 9, height = 7, dpi = 300)
plots$CytoTRACE2_Relative_UMAP
ggsave("CytoTRACE2_Relative_UMAP.pdf",  width = 9, height = 7, dpi = 300)
plots$Phenotype_UMAP
ggsave("Phenotype_UMAP.pdf",  width = 9, height = 7, dpi = 300)
plots$CytoTRACE2_Boxplot_byPheno
ggsave("CytoTRACE2_Boxplot_byPheno.pdf",  width = 5, height = 3, dpi = 300)


library(paletteer)
library(Seurat)
library(monocle3)
library(dplyr)
library(BiocParallel)
register(MulticoreParam(workers = 8, progressbar = TRUE))
scRNA <- sub_data
Idents(scRNA) <- scRNA$cell_type
levels(Idents(scRNA)) 
scRNA$celltype <- factor(scRNA$cell_type,levels = c("Dendritic cells","Macrophages","Microglial cells", "Monocytes","Neutrophils" ))
Idents(scRNA) <- scRNA$celltype

## 提取数据
expression_matrix <- GetAssayData(scRNA, assay = 'RNA',slot = 'counts')
cell_metadata <- scRNA@meta.data 
gene_annotation <- data.frame(gene_short_name = rownames(expression_matrix))
rownames(gene_annotation) <- rownames(expression_matrix)
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
# 归一化/预处理数据 这一步使用的PCA分析，dim数代表纳入的PCA数量
cds <- preprocess_cds(cds, num_dim = 25)
# 这个函数用于确认设定的dim数是否足够代表主要变异
plot_pc_variance_explained(cds)

# 可选(去批次处理)
#cds <- align_cds(cds, num_dim = 100, alignment_group = "GSE_num")
# 降维聚类，可选择UMAP、PCA或者TSNE
cds <- reduce_dimension(cds,reduction_method='UMAP',preprocess_method = 'PCA')

plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "celltype",group_label_size = 5,rasterize = FALSE)
#基因/轨迹共定位  这里可以看基因和轨迹的共定位情况
modelgenes <- "Spp1"
plot_cells(cds,
           cell_size=1.5,group_label_size=6,
           genes=modelgenes,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)
ggsave("genes_show.pdf",width = 12,height = 12)

#6.细胞聚类
cds <- cluster_cells(cds) #cluster your cells
plot_cells(cds, color_cells_by = "celltype",group_label_size = 3.5,rasterize = FALSE)
#7.轨迹分析# 轨迹推断
cds <- learn_graph(cds,verbose=T,
                   use_partition=T, #默认是T，T时是顾及全局的情况
)
plot_cells(cds, 
           color_cells_by = 'celltype',
           label_groups_by_cluster=FALSE,
           cell_size=1,group_label_size=4,
           label_leaves= T, # 是否显示不同细胞结局
           label_branch_points=T, # 是否显示不同的分支节点
           trajectory_graph_color='#023858',
           trajectory_graph_segment_size = 1)
#8.定义起点-轨迹可视化
#定义root cell, 推断拟时方向
# 网页自定
cds <- order_cells(cds)  
# 代码定
# a helper function to identify the root principal points:
## 指定初始细胞群(ZSCAN12+T)和细胞亚群的列名(celltype)
get_earliest_principal_node <- function(cds, time_bin="Dendritic cells"){
  cell_ids <- which(colData(cds)[, "celltype"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
## 2D 可视化
plot_cells(cds, label_cell_groups = T, 
           color_cells_by = "pseudotime", 
           label_branch_points = F, 
           label_roots =F, # 显示根
           label_leaves =F,
           graph_label_size = 0, 
           cell_size=2, 
           trajectory_graph_color='black',
           trajectory_graph_segment_size = 2)




#9.轨迹差异基因分析
# 提取不同轨迹的差异基因/并选择前12个
# neighbor_graph="principal_graph"提取轨迹上相似位置是否有相关的表达
trace_genes <- graph_test(cds, 
                          neighbor_graph = "principal_graph", 
                          cores = 10)
saveRDS(trace_genes,file = 'trace_genes.rds')
trace_genes<-read_rds('trace_genes.rds')
track_genes_sig <- trace_genes %>%
  top_n(n=12, morans_I) %>% 
  pull(gene_short_name) %>% 
  as.character()

# 差异基因绘制######################################
levels(Idents(scRNA)) #打出来细胞类型供复制
# [1] "ZSCAN12+T" "Naive T"   "Th1"       "Tm"        "ZNF793+T"  "ELK4+T"    "Th17"     
# [8] "Treg"
lineage_cds <- cds[rowData(cds)$gene_short_name %in% track_genes_sig,
                   colData(cds)$celltype %in% c("Dendritic cells","Macrophages","Microglial cells", "Monocytes","Neutrophils" )]
#lineage_cds <- order_cells(lineage_cds)
plot_genes_in_pseudotime(lineage_cds,
                         color_cells_by="celltype",
                         min_expr=0.5)
ggsave("degs_monocles.pdf",width = 8,height = 14)

# 细胞映射
plot_cells(cds, genes= track_genes_sig,
           cell_size=1, 
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)
ggsave("track_genes_sig.pdf",width = 12,height = 9)


library(ClusterGVis)
genes <- row.names(subset(trace_genes, q_value < 0.001 & morans_I > 0.35))
mat <- pre_pseudotime_matrix(cds_obj = cds,gene_list = genes)
mat <- as.data.frame(mat)
head(mat[1:5,1:5])

# kmeans聚类
ck <- ClusterGVis::clusterData(obj = mat,
                               cluster.method = "kmeans",
                               cluster.num = 4)

# 富集分析
library(org.Hs.eg.db)##人是org.Hs.eg.db
library(org.Mm.eg.db)
enrich <- enrichCluster(object = ck,
                        OrgDb = org.Mm.eg.db, ##人是org.Hs.eg.db
                        type = "BP", # 可以自己定义
                        pvalueCutoff = 0.01,
                        topn = 5, # topn设置为NULL可获得全部结果
                        seed = 123)
options(bitmapType = "cairo")  # 使用Cairo图形设备
memory.limit(size = 160000000)     # 增加内存至16GB（根据系统调整）

# 添加行注释
#pdf('monocle3.pdf',height = 10,width = 12,onefile = F)
visCluster(object = ck,
           plot.type = "both", 
           add.sampleanno = F,
           markGenes.side = "left",
           markGenes = sample(rownames(mat),30,replace = F),
           ht.col.list = list(col_range = c(-2,-1.5,-1,-0.5,0,0.5,1,1.5,2),
                              col_color = c("#0D0887FF","#42049EFF","#6A00A8FF",
                                            "#900DA4FF","white","#E16462FF",
                                            "#FCA636FF","#FCCE25FF","#F0F921FF")),
           genes.gp = c('italic',fontsize = 12,col = "black"),
           line.side = "left",
           annoTerm.data = enrich
           )
dev.off()

