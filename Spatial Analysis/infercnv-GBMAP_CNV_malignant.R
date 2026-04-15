# 加载必要的包
library(infercnv)
library(Seurat)
library(dplyr)
library(tidyverse)
library(RColorBrewer) # 提供Brewer调色板
library(ggsci)
mycol <- c(pal_d3()(7),pal_aaas()(7),pal_uchicago()(7),pal_jama()(7))

# 读取参考数据集
scRef <- readRDS("GBMAP_ssgsea.RDS")
DimPlot(scRef, group.by  = "iCNV",raster=FALSE,cols =mycol)

# 创建细胞注释数据框
cell_annotations <- data.frame(
  cell_id = colnames(scRef),
  cell_type = scRef@meta.data$annotation_level_4
)
getwd()
setwd("D:/ProgramData/Glioma-immunelncRNA/spatial-intergrate-gbmap/")
# 定义参考细胞类型（非恶性细胞）
# 基于您的细胞类型注释，推荐以下细胞作为reference：
reference_celltypes <- c("B_cell", "Plasma_B", "CD4_rest", "CD4_INF", 
                         "CD8_cytotoxic", "CD8_EM", "CD8_NK_sig", "NK",
                         "cDC1", "cDC2", "DC1", "DC2", "DC3", "pDC", "Mast",
                         "Astrocyte", "Oligodendrocyte", "Pericyte", 
                         "Endo_arterial", "Endo_capilar", "SMC", "VLMC")

# 提取恶性细胞（观察组）- 基于常见的胶质瘤恶性表型
malignant_celltypes <- c("AC_like", "AC_like_Prolif", "MES_like_hypoxia_independent",
                         "MES_like_hypoxia_MHC", "NPC_like_neural", "NPC_like_OPC",
                         "NPC_like_Prolif", "OPC_like", "OPC_like_Prolif")

# 筛选需要分析的细胞（参考细胞 + 恶性细胞）
selected_cells <- cell_annotations %>%
  filter(cell_type %in% c(reference_celltypes, malignant_celltypes))

# 创建用于inferCNV的细胞注释文件
infercnv_annotations <- data.frame(
  cell_id = selected_cells$cell_id,
  cell_type = selected_cells$cell_type,
  stringsAsFactors = FALSE
)

# 保存注释文件（临时）

write.table(infercnv_annotations, "./GBMAP_CNV_malignant/groupFiles.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# 提取原始计数矩阵（只包含选择的细胞）
exprMatrix <- GetAssayData(scRef, slot = "counts")[, selected_cells$cell_id]

# 创建inferCNV对象

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=exprMatrix,#单细胞Counts矩阵 
                                    annotations_file="./GBMAP_CNV_malignant/groupFiles.txt",#分组文件(细胞) 
                                    delim="\t", 
                                    gene_order_file= "./GBMAP_CNV_malignant/cnv_ref/hg38_gencode_v27.txt",#基因位置文件 
                                    ref_group_names= reference_celltypes #这里使用Control 
)


# 运行inferCNV分析
infercnv_obj = infercnv::run(infercnv_obj,#inferCNV对象 
                             cutoff=0.1, # use 1 for smart-seq, 0.1 for 10x-genomics 
                             out_dir= './GBMAP_CNV_malignant/' , # dir is auto-created for storing outputs 
                             cluster_by_groups= T, # cluster_by_groups：先区分细胞来源，再做层次聚类 
                             hclust_method="ward.D2", 
                             plot_steps=F, analysis_mode = "subclusters",# 分析模式
                             num_threads = 4, 
                             write_expr_matrix = T, #重要！！ # 输出表达矩阵
                             output_format = "pdf",  # 使用 PDF 而非默认的 PNG（更节省内存）
                             res = 200, useRaster = TRUE
)
infercnv_obj <- readRDS("GBMAP_CNV_malignant/run.final.infercnv_obj") # 请将路径替换为你的实际输出目录
str(infercnv_obj) # 简单查看对象结构
# 3. 克隆识别 ----
# 获取CNV矩阵
cnv_matrix <- infercnv_obj@expr.data
setwd('./GBMAP_CNV_malignant/')
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

infercnv::plot_cnv(
  infercnv_obj,
  out_dir = "GBMAP_CNV_malignant",
  title = "inferCNV Analysis - GBM Malignant Cells",
  obs_title = "Putative Malignant Cells",
  ref_title = "Reference Cells (Non-malignant)",
  color_scheme = rev(brewer.pal(9, "RdBu")),  # 红蓝色系
  xlab = "Genomic Position",
  ylab = "Cells",
  cluster_by_groups = TRUE
)
