.libPaths(c('~/R/x86_64-pc-linux-gnu-library/4.3/'))
setwd("/home/data/t060202/Spatial-GSVA-duijian/Spatial_Glioma-CELL/intergreted/")
# 加载所需包
library(Seurat)
library(patchwork)
library(ggplot2)
library(bigmemory)
library(doMC)
library(patchwork)
library(future)
library(presto)
library(stringr)
library(harmony)
library(dplyr)
library(HGNChelper)
library(future)
library(ggsci)
library(quadprog)
library(singlet)
##1读取数据
plan("multisession", workers = 16)
plan()
#设置可用的内存
options(future.globals.maxSize = 100* 1024^3)
###### Per sample Leiden clustering##############################
mycol <- c(pal_d3()(7),pal_aaas()(7),pal_uchicago()(7),pal_jama()(7))

samples_name <- (read.delim("GBM_samples.txt", header = FALSE, sep = "\t"))$V1
# 设置参数
metadata <- read.csv("visium_metadata.csv", header = TRUE)          # 替换为你的元数据文件路径
# 合并前两列为新列名（示例：SampleID和Batch列）
new_row_names <- paste(metadata$sample, metadata$spot_id, sep = "_")
length(new_row_names)
length(rownames(metadata))
# 更新列名（保留原数据内容）
rownames(metadata) <- new_row_names
metadata <- metadata[, -c(1)]  

spatial.list = list()
for(i in 1:length(samples_name)){
  spatial <- Load10X_Spatial(data.dir = paste("../general/GBM_data/",samples_name[i],"/outs", sep = ""),
                             slice = samples_name[i],
                             assay = "Spatial",
                             filter.matrix = TRUE,
                             to.upper = FALSE)
  spatial$orig.ident <- samples_name[i]
  spatial <- RenameCells(spatial, add.cell.id = samples_name[i])
  spatial [["percent.mt"]] <- PercentageFeatureSet(spatial, pattern = "^MT[-]")
  spatial.list[[samples_name[i]]] <- spatial
}
# 2. 标准化和变量基因选择（各样本独立处理）
seurat_list <- lapply(spatial.list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, nfeatures = 2000)
})

# 3. 数据整合
# 选择整合用基因
features <- SelectIntegrationFeatures(seurat_list)

# 找整合锚点
anchors <- FindIntegrationAnchors(
  object.list = seurat_list,
  anchor.features = features
)

# 创建整合数据集
integrated <- IntegrateData(anchors)

# 4. 处理负值问题（NMF需要非负输入）
# 获取整合后的数据矩阵
integrated_data <- GetAssayData(integrated, assay = "integrated", slot = "data")

# 将负值设为0（根据需求调整）
integrated_data[integrated_data < 0] <- 0
integrated <- SetAssayData(integrated, 
                           assay = "integrated",
                           slot = "data",
                           new.data = integrated_data)

# 5. 运行NMF
set.seed(42)
integrated <- RunNMF(
  object = integrated,
  assay = "integrated",
  n.cores = 4,      # 并行计算
  tol = 1e-5,       # 收敛阈值
  maxit = 100       # 最大迭代次数
)

# 6. 可视化结果
# 显示交叉验证结果
RankPlot(integrated)

# 可视化前6个因子
selected_factors <- paste0("NMF_", 1:6)
MapFeatures(integrated, 
            features = selected_factors,
            colors = viridis::magma(n = 11, direction = -1),
            pt_size = 1.5)

# 查看基因贡献度
PlotFeatureLoadings(integrated,
                    dims = 1:3,
                    nfeatures = 20,
                    mode = "heatmap")