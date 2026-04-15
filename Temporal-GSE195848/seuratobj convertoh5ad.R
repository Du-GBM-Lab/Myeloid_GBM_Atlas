library(scCustomize)
library(Seurat)
library(reticulate)
reticulate::use_condaenv("sctour", required = TRUE)
use_condaenv('sctour')

meyloid <- readRDS("meyloid_SeuratObj.RDS")
Assays(meyloid)
DefaultAssay(meyloid) <- 'RNA'
DimPlot(meyloid,label = T)
# 转化为h5ad文件
as.anndata(x = meyloid, file_path = "./", 
           file_name = "meyloid.h5ad", assay = "RNA", main_layer = "counts", 
           other_layers = c("data"), transfer_dimreduc = TRUE, verbose = TRUE)
