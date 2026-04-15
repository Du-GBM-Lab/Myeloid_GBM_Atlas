library(Seurat)
library(SeuratData)
library(SeuratObject)
library(harmony)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(forcats)
library(SingleCellExperiment)
library(RColorBrewer)
library(SeuratData)
library(RcppML)
library(dittoSeq)
library(future)
library(randomcoloR)
library(ggsci)
library(scCustomize)
library(infercnv) 
library(miloR)          # 邻域分析
library(pheatmap)
library(scran)
library(scater)
library(patchwork)
mycol <- c(pal_d3()(7),pal_aaas()(7),pal_uchicago()(7),pal_jama()(7))
setwd('D:/ProgramData/Glioma-immunelncRNA/Temporal-GSE195848/')
meyloid_SeuratObj<-readRDS('meyloid_SeuratObj.RDS')
table(meyloid_SeuratObj$cell_type)
table(meyloid_SeuratObj$stage)
seurat_obj = meyloid_SeuratObj[,meyloid_SeuratObj$stage %in% c("Normal","Medium")]
dim(seurat_obj)
table(seurat_obj$stage)

DimPlot(seurat_obj,split.by = 'stage',group.by = 'cell_type')
exprMatrix <- as.matrix(GetAssayData(seurat_obj, slot = 'counts'))
# 提取细胞类型信息 
cellAnnota <- subset(seurat_obj@meta.data, select = c('stage')) 
write.table(cellAnnota, "./CNV_Medium/groupFiles.txt", sep = '\t', col.names = FALSE)


#---------------------创建InferCNV对象------------------------- 
#ref_group_names参数根据细胞注释文件填写，在示例中，这两种细胞是非恶性细胞，所以作为参照； 
#ref_group_names=NULL，则会选用所有细胞的平均表达作为参照，这种模式下，最好能确保细胞内在有足够的差异 
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=exprMatrix,#单细胞Counts矩阵 
                                    annotations_file="./CNV_Medium//groupFiles.txt",#分组文件(细胞) 
                                    delim="\t", 
                                    gene_order_file= "./CNV_Medium/cnv_ref/mouse_gencode.GRCm39.vM32.basic.annotation.by_gene_name.infercnv_positions",#基因位置文件 
                                    ref_group_names= c('Normal') #这里使用Control 
)
gc() 


# 检查内存限制
memory.limit()
rm(exprMatrix,seurat_obj)
infercnv_obj = infercnv::run(infercnv_obj,#inferCNV对象 
                             cutoff=0.1, # use 1 for smart-seq, 0.1 for 10x-genomics 
                             out_dir= './CNV_Medium/' , # dir is auto-created for storing outputs 
                             cluster_by_groups= T, # cluster_by_groups：先区分细胞来源，再做层次聚类 
                             hclust_method="ward.D2", 
                             plot_steps=F, 
                             num_threads = 30, 
                             write_expr_matrix = T, #重要！！ 
                             output_format = "pdf",  # 使用 PDF 而非默认的 PNG（更节省内存）
                             res = 200, useRaster = TRUE
)
# 3. 克隆识别 ----
# 获取CNV矩阵
cnv_matrix <- infercnv_obj@expr.data
setwd('./CNV_Medium/')
getwd()
load('infercnv_obj-hclust.RData')
# 层次聚类定义克隆
dist_matrix <- dist(t(cnv_matrix))
hclust_results <- hclust(dist_matrix, method = "ward.D2")
save(infercnv_obj,hclust_results, file = 'infercnv_obj-hclust.RData')

clone_ids <- cutree(hclust_results, k = 5)  # 根据轮廓系数调整k值

# 添加到细胞注释
meyloid_SeuratObj$clone_id <- clone_ids[colnames(meyloid_SeuratObj)]
table(meyloid_SeuratObj$clone_id)
# 4. 邻域分析（miloR） ----
##转换对象
meyloid_sce <- as.SingleCellExperiment(meyloid_SeuratObj)
# 创建Milo对象
milo_obj <- Milo(meyloid_sce)

reducedDim(milo_obj, "UMAP") <- reducedDim(meyloid_sce, "UMAP")
plotUMAP(milo_obj, group_by="cell_type") + plotUMAP(milo_obj, colour_by="stage")
# 示例代码
plotUMAP(milo_obj, colour_by="cell_type") + 
  scale_color_npg() +  # Nature 风格配色
  theme_minimal()

plotUMAP(milo_obj, colour_by="stage") + 
  scale_color_brewer(palette = "Set3")  # RColorBrewer 调色板
# 构建邻域图
milo_obj <- buildGraph(milo_obj, k = 20, d = 10)  # k=邻域大小, d=降维维度

# 定义邻域
milo_obj <- makeNhoods(milo_obj, prop = 0.1)  # prop=采样比例
milo_obj <- countCells(milo_obj, meta.data = as.data.frame(colData(milo_obj)), sample="sample")

# 计算邻域统计量
milo_obj <- calcNhoodDistance(milo_obj, d = 10)
#3.5定义实验设计
milo_design <- data.frame(colData(milo_obj))[,c("sample", "stage", "batch")]

## 将批次信息从整数转换为因子
milo_design$batch <- as.factor(milo_design$batch) 
milo_design <- distinct(milo_design)
rownames(milo_design) <- milo_design$sample

milo_design

#3.7 测试
da_results <- testNhoods(milo_obj, design = ~ batch + stage, design.df = milo_design)
head(da_results)
ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)
ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1) ## 标记显著性阈值（10% FDR）

milo_obj <- buildNhoodGraph(milo_obj)
reducedDimNames(milo_obj)
## 绘制单细胞UMAP
umap_pl <- plotReducedDim(milo_obj, dimred = "UMAP", colour_by="stage", text_by = "cell_type", 
                          text_size = 3, point_size=0.5) +
  guides(fill="none")

## 绘制邻域图
nh_graph_pl <- plotNhoodGraphDA(milo_obj, da_results, layout="UMAP",alpha=0.1) 

umap_pl + nh_graph_pl + plot_layout(guides="collect")
ggsave('patchwork-miloR-Medium.pdf',width = 33.6, height = 12.8, units = "cm")
ggsave('patchwork-miloR-Medium1.pdf',width = 14.4, height = 7.4, units = "in")



# 5. 可视化 ----
# Fig.3a左：克隆热图
# 准备热图数据
heatmap_data <- table(meyloid_SeuratObj$clone_id, 
                      meyloid_SeuratObj$cell_type)  # cell_state需替换为你的实际注释列名
table(meyloid_SeuratObj$cell_type)
pheatmap<-pheatmap(heatmap_data,
                   color = colorRampPalette(c("white", "firebrick"))(50),
                   cluster_rows = TRUE,
                   cluster_cols = TRUE,
                   show_rownames = TRUE,angle_col = '45',
                   main = "Clone Composition by Cell Type")
ggsave('Clone Composition by Cell_Medium.pdf',plot = pheatmap,width = 6, height = 4, units = "in")

##各个细胞亚群内的邻域丰度分布
da_results <- annotateNhoods(milo_obj, da_results, coldata_col = "cell_type")
da_results$cell_type <- factor(da_results$cell_type, 
                              levels=c('Macrophages', 'Microglial cells', 'Monocytes', 'Neutrophils', 'Dendritic cells'), 
                              ordered=T)

write.csv(da_results, 'result_Medium.csv', quote=F)

da_plot <- da_results[da_results$SpatialFDR < 0.1, ]
p<-plotDAbeeswarm(da_plot, group.by = "cell_type") + 
  scale_color_gradient2(high='red', mid='white', low='blue') + xlab('')
ggsave('plotDAbeeswarm_medium.pdf',plot = p,width = 6, height = 4, units = "in")

