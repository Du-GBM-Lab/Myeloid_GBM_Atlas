.libPaths()
.libPaths(c('~/R/x86_64-pc-linux-gnu-library/4.3/'))
.libPaths(c('~/R/x86_64-pc-linux-gnu-library/4.4/'))
getwd()
setwd("/home/data/t060202/Spatial-GSVA-duijian/GBMAP-data/")
getwd()
library(sceasy)
library(reticulate)
library(Seurat)
library(bigmemory)
library(future)
library(HGNChelper)
library(future)
library(quadprog)
library(schard)
library(BPCells)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(harmony)
library(plot1cell)
library(ggsci)
library(quadprog)
library(RColorBrewer)
library(scCustomize)
library(Azimuth)
library(SeuratDisk)
library(qs)
library(harmony)
library(scCustomize)
library(infercnv) 
library(miloR)          # 邻域分析
library(pheatmap)
plan("multisession", workers = 4)
plan()
#设置可用的内存
options(future.globals.maxSize = 1000000* 1024^5)
###### Per sample Leiden clustering##############################
mycol <- c(pal_d3()(7),pal_aaas()(7),pal_uchicago()(7),pal_jama()(7),pal_gsea()(7),pal_cosmic()(7),pal_flatui()(7),pal_igv()(7))
obj <- qread("obj.qs")

meyloid_SeuratObj1<-obj
meyloid_SeuratObj1[["percent.mt"]] <- PercentageFeatureSet(meyloid_SeuratObj1, pattern = "^MT")
meyloid_SeuratObj1 <- subset(meyloid_SeuratObj1,subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & percent.mt < 15)  # 调整阈值
VlnPlot(meyloid_SeuratObj1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
DefaultAssay(meyloid_SeuratObj1) <- "RNA"
meyloid_SeuratObj1 <- NormalizeData(meyloid_SeuratObj1,normalization.method = "LogNormalize",scale.factor = 1e4,verbose = FALSE)
meyloid_SeuratObj1[["RNA"]]@meta.features <- data.frame(row.names = rownames(meyloid_SeuratObj1[["RNA"]]))
meyloid_SeuratObj1 <- FindVariableFeatures(meyloid_SeuratObj1, selection.method = "vst",nfeatures = 3000)
meyloid_SeuratObj1 <- ScaleData(meyloid_SeuratObj1,vars.to.regress = c("nFeature_RNA", "percent.mt"))
meyloid_SeuratObj1 <- RunPCA(meyloid_SeuratObj1, features = VariableFeatures(object = meyloid_SeuratObj1))
meyloid_SeuratObj1 <- FindNeighbors(meyloid_SeuratObj1, dims = 1:20)
meyloid_SeuratObj1 <- FindClusters(meyloid_SeuratObj1,resolution = seq(0.1, 1, by = 0.1), 
                                  verbose = FALSE)  # 分辨率按需调整
library(clustree)
library(patchwork)
p1 <- clustree(meyloid_SeuratObj1, prefix = 'RNA_snn_res.') + coord_flip();p1
#这里的RNA_snn_res后面的数值是可以修改的
p2 <- DimPlot(meyloid_SeuratObj1, group.by = 'RNA_snn_res.0.8', label = T);p2 
p2 <- DimPlot(meyloid_SeuratObj1, group.by = 'RNA_snn_res.0.5', label = T);p2 
Tree_1 <- p1 + p2 + plot_layout(widths = c(3, 1));print(Tree_1)
#UMAP/tSNE可视化前先确定一个想要的reslution值
#这里重新跑一遍之后后面就会按照新跑的reslution值进行分析
meyloid_SeuratObj1 <- FindClusters(meyloid_SeuratObj1, resolution = 0.5, verbose = FALSE)
meyloid_SeuratObj1 <- RunUMAP(meyloid_SeuratObj1, dims = 1:20)
meyloid_SeuratObj1 <- RunTSNE(meyloid_SeuratObj1, dims = 1:20)

DimPlot(meyloid_SeuratObj1, reduction = "umap", label = TRUE,group.by = 'annotation_level_3')
DimPlot(meyloid_SeuratObj1, reduction = "tsne", label = TRUE,cols = mycol,group.by = 'annotation_level_3')

qsave(meyloid_SeuratObj1, "meyloid_exGBmap-stand.qs")
meyloid_SeuratObj1 <- qread("meyloid_exGBmap-stand.qs")
meyloid_SeuratObj1@reductions$Xumap_ <- NULL
meyloid_SeuratObj1@meta.data[["gbmap"]] <- NULL
meyloid_SeuratObj1@meta.data[["observation_joinid"]] <- NULL

colnames(meyloid_SeuratObj1@meta.data)
DimPlot(meyloid_SeuratObj1,reduction = 'Xumap_',group.by = 'annotation_level_3')
DimPlot(meyloid_SeuratObj1,reduction = 'pca',group.by = 'annotation_level_3')
DimPlot(meyloid_SeuratObj1,reduction = 'umap',group.by = 'annotation_level_3')
DimPlot(meyloid_SeuratObj1,reduction = 'tsne',group.by = 'disease')
DimPlot(obj,reduction = 'Xumap_',group.by = 'observation_joinid')
DimPlot(obj,reduction = 'Xumap_',group.by = 'annotation_level_3')
DimPlot(obj,reduction = 'Xumap_',group.by = 'annotation_level_3')

qsave(meyloid_SeuratObj1, "meyloid-stand_upload.qs")
