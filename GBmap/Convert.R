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
use_condaenv('sceasy')
loompy <- reticulate::import('loompy')
use_condaenv('leiden')
myeloid_subset@version  # 查看对象版本
packageVersion("Seurat") # 查看当前安装的Seurat版本
sceasy::convertFormat("myeloid_subset.h5ad",from="anndata",to="seurat",outFile='myeloid_subset.rds')

myeloid_subset = schard::h5ad2seurat('myeloid_subset.h5ad',use.raw=T)
colnames(myeloid_subset@meta.data)
DimPlot(myeloid_subset,group.by = "annotation_level_3")
head(rownames(myeloid_subset))
myeloid_subset_v5 <- UpdateSeuratObject(myeloid_subset)
# 升级对象后，显式创建 Assay5
myeloid_subset_v5[["RNA5"]] <- as(object = myeloid_subset_v5[["RNA"]], Class = "Assay5")
myeloid_subset_v5[["RNA"]] <- NULL
# 设置新创建的 Assay5 为默认分析对象
DefaultAssay(myeloid_subset_v5) <- "RNA5"
myeloid.Azimuth <- Azimuth:::ConvertEnsembleToSymbol(mat = myeloid_subset_v5, species = "human")
myeloid.Azimuth[["RNA"]] <- NULL
head(rownames(myeloid.Azimuth))
library(org.Hs.eg.db)
head(rownames(myeloid_subset))
ids=select(org.Hs.eg.db,keys = rownames(myeloid_subset),
           columns = c('ENSEMBL','SYMBOL'),
           keytype = 'ENSEMBL')
head(ids)
dim(ids) # [1] 16428 
ids=na.omit(ids)
dim(ids) # [1] 15504 
length(unique(ids$SYMBOL)) # [1] 15494 
# 这里的关系超级乱，互相之间都不是一对一 
# 凡是混乱的ID一律删除即可
ids=ids[!duplicated(ids$SYMBOL),]
ids=ids[!duplicated(ids$ENSEMBL),]
pos=match(ids$ENSEMBL,rownames(myeloid_subset) )
myeloid_subset=myeloid_subset[pos,]
myeloid_subset # 15393 features
# https://github.com/satijalab/seurat/issues/2617
# RenameGenesSeurat  ------------------------------------------------------------------------------------
RenameGenesSeurat <- function(obj , 
                              newnames ) { 
  # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. 
  # It only changes obj@assays$RNA@counts, @data and @scale.data.
  print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
  RNA <- obj@assays$RNA
  
  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
    if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays$RNA <- RNA
  return(obj)
}
obj=RenameGenesSeurat(obj = myeloid_subset, 
                      newnames = ids$SYMBOL)

myeloid_subset
obj
myeloid.Azimuth
head(rownames(obj))
head(obj@assays$RNA@counts)
DimPlot(obj,group.by = 'annotation_level_2',cols = mycol,raster=FALSE)
table(obj$annotation_level_2)
table(obj$annotation_level_4)
qsave(obj, "obj.qs")
obj <- qread("obj.qs")



# 提取两个对象的基因名
genes_myeloid <- rownames(myeloid.Azimuth)
genes_obj <- rownames(obj)

# 检查基因数量差异
cat("myeloid.Azimuth 基因数:", length(genes_myeloid), "\n")
cat("obj 基因数:", length(genes_obj), "\n")

unique_to_myeloid <- setdiff(genes_myeloid, genes_obj)
cat("仅存在于 myeloid.Azimuth 的基因数:", length(unique_to_myeloid), "\n")
head(unique_to_myeloid)  # 预览前几个基因

# 查找仅存在于 obj 的基因
unique_to_obj <- setdiff(genes_obj, genes_myeloid)
cat("仅存在于 obj 的基因数:", length(unique_to_obj), "\n")
head(unique_to_obj)

# 检查共有基因数
common_genes <- intersect(genes_myeloid, genes_obj)
cat("共有基因数:", length(common_genes), "\n")

meyloid_SeuratObj<-obj
meyloid_SeuratObj[["percent.mt"]] <- PercentageFeatureSet(meyloid_SeuratObj, pattern = "^MT")
meyloid_SeuratObj <- subset(meyloid_SeuratObj,subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & percent.mt < 15)  # 调整阈值
VlnPlot(meyloid_SeuratObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
DefaultAssay(meyloid_SeuratObj) <- "RNA"
meyloid_SeuratObj <- NormalizeData(meyloid_SeuratObj,normalization.method = "LogNormalize",scale.factor = 1e4,verbose = FALSE)
meyloid_SeuratObj[["RNA"]]@meta.features <- data.frame(row.names = rownames(meyloid_SeuratObj[["RNA"]]))
meyloid_SeuratObj <- FindVariableFeatures(meyloid_SeuratObj, selection.method = "vst",nfeatures = 3000)
meyloid_SeuratObj <- ScaleData(meyloid_SeuratObj,vars.to.regress = c("nFeature_RNA", "percent.mt"))
meyloid_SeuratObj <- RunPCA(meyloid_SeuratObj, features = VariableFeatures(object = meyloid_SeuratObj))
meyloid_SeuratObj <- RunHarmony(object = meyloid_SeuratObj,group.by.vars = c('author'), dims.use = 1:20, max.iter.harmony = 50)
ElbowPlot(meyloid_SeuratObj)  # 选择主成分数（如dims=1:15）
names(seuratObj@reductions)
meyloid_SeuratObj <- RunUMAP(meyloid_SeuratObj, dims = 1:15, reduction = "harmony")
meyloid_SeuratObj <- RunTSNE(meyloid_SeuratObj, dims = 1:15, reduction = "harmony")

DimPlot(seuratObj, reduction = "umap", label = TRUE) 






meyloid_SeuratObj <- FindNeighbors(meyloid_SeuratObj, dims = 1:20)
meyloid_SeuratObj <- FindClusters(meyloid_SeuratObj,resolution = seq(0.1, 1, by = 0.1), 
                                  verbose = FALSE)  # 分辨率按需调整
cluster_harmony <- DimPlot(meyloid_SeuratObj, reduction = "harmony");print(cluster_harmony) #单样本需要pca
library(clustree)
library(patchwork)
p1 <- clustree(meyloid_SeuratObj, prefix = 'RNA_snn_res.') + coord_flip();p1
#这里的RNA_snn_res后面的数值是可以修改的
p2 <- DimPlot(meyloid_SeuratObj, group.by = 'RNA_snn_res.0.8', label = T);p2 
Tree_1 <- p1 + p2 + plot_layout(widths = c(3, 1));print(Tree_1)
#UMAP/tSNE可视化前先确定一个想要的reslution值
#这里重新跑一遍之后后面就会按照新跑的reslution值进行分析
meyloid_SeuratObj <- FindClusters(meyloid_SeuratObj, resolution = 0.8, verbose = FALSE)
meyloid_SeuratObj <- RunUMAP(meyloid_SeuratObj, dims = 1:20)
meyloid_SeuratObj <- JoinLayers(meyloid_SeuratObj)
DimPlot(meyloid_SeuratObj, reduction = "umap", label = TRUE,cols = mycol,group.by = 'cell_type')
DimPlot(meyloid_SeuratObj, reduction = "umap", label = TRUE,cols = mycol,group.by = 'cell_type',split.by = 'stage')
saveRDS(meyloid_SeuratObj,file = 'meyloid_SeuratObj.RDS')
meyloid_SeuratObj<-readRDS('meyloid_SeuratObj.RDS')












