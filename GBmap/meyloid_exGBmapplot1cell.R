.libPaths()
.libPaths(c('~/R/x86_64-pc-linux-gnu-library/4.3/'))
getwd()
setwd("/home/data/t060202/Spatial-GSVA-duijian/GBMAP-data/")
library(Seurat)
library(monocle)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(knitr)
lapply(c("dplyr", "HGNChelper", "openxlsx"), library, character.only = T)
library(qs)
library(ggsci)
library(plot1cell)
library(ggsci)
library(quadprog)
library(RColorBrewer)
library(scCustomize)
library(presto)
plan("multisession", workers = 16)
plan()
#设置可用的内存
options(future.globals.maxSize = 10000* 1024^5)
###### Per sample Leiden clustering##############################
mycol <- c(pal_d3()(7),pal_aaas()(7),pal_uchicago()(7),pal_jama()(7),pal_gsea()(7),pal_cosmic()(7),pal_flatui()(7),pal_igv()(7))

meyloid_SeuratObj<-qread('meyloid_exGBmap-stand.qs')

DimPlot(meyloid_SeuratObj, reduction = "tsne", label = TRUE,cols = mycol,group.by = 'annotation_level_3',raster=FALSE)
DimPlot(meyloid_SeuratObj, reduction = "tsne", label = TRUE,cols = mycol,group.by = 'annotation_level_2',raster=FALSE)
DimPlot(TAM_MG_SeuratObj, reduction = "tsne", label = TRUE,cols = mycol,group.by = 'annotation_level_4',raster=FALSE,, label.size = 6,repel = T)
DimPlot(meyloid_SeuratObj,group.by = 'annotation_level_3',cols = mycol,raster=FALSE,, label.size = 6,repel = T,reduction = "umap")
get_metadata <- function(
    seu_obj, 
    reductions = "tsne", 
    coord_scale = 0.8, 
    color
){
  metadata<-seu_obj@meta.data
  metadata$Cluster<-seu_obj@active.ident
  metadata$dim1<-as.numeric(seu_obj[[reductions]]@cell.embeddings[,1])
  metadata$dim2<-as.numeric(seu_obj[[reductions]]@cell.embeddings[,2])
  metadata$x<-transform_coordinates(metadata$dim1, zoom = coord_scale)
  metadata$y<-transform_coordinates(metadata$dim2, zoom = coord_scale)
  color_df<-data.frame(Cluster=levels(seu_obj), Colors=color)
  cellnames<-rownames(metadata)
  metadata$cells<-rownames(metadata)
  metadata<-merge(metadata, color_df, by='Cluster')
  rownames(metadata)<-metadata$cells
  metadata<-metadata[cellnames,]
  metadata
}

prepare_circlize_data <- function(
    seu_obj, 
    scale =0.8
){
  celltypes<-levels(seu_obj)
  cell_colors <- scales::hue_pal()(length(celltypes))
  data_plot <- get_metadata(seu_obj, color = cell_colors, coord_scale = scale)
  data_plot <- cell_order(data_plot)
  data_plot$x_polar2 <- log10(data_plot$x_polar)
  data_plot
}


###Prepare data for ploting

meyloid_SeuratObj$Cluster <- meyloid_SeuratObj$annotation_level_3
meyloid_SeuratObj <- SetIdent(meyloid_SeuratObj, value = "Cluster")
###Prepare data for ploting
circ_data <- prepare_circlize_data(meyloid_SeuratObj, scale = 0.75)
set.seed(1234)
levels(meyloid_SeuratObj$Cluster)
circ_data$Cluster <- as.factor(circ_data$Cluster)
levels(meyloid_SeuratObj$annotation_level_3)
levels(circ_data$Cluster)
colnames(circ_data)
circ_data$Cluster <- factor(x = circ_data$Cluster, 
                            levels = c('DC','Mast','Mono','TAM-BDM','TAM-MG','Neutrophil'))
col = c(
  "TAM - BDM" = "#FF6B6B",  # 亮红色替代紫色
  "TAM - MG" = "#4CC9F0",   # 亮蓝色
  "Mono" = "#8AC926",       # 苹果绿
  "DC" = "#7209B7",         # 深紫色
  "Mast" = "#F72585",        # 玫粉色
  "Neutrophil"= "#FFA726" 
)
tissue_colors<-rand_color(length(names(table(meyloid_SeuratObj$tissue))))
TERT_colors<-rand_color(length(names(table(meyloid_SeuratObj$TERT))))

location_colors<-rand_color(length(names(table(meyloid_SeuratObj$location))))



library(circlize)
png(filename =  'circlize_plot.png', width = 8, height = 8, units = 'in', res = 300)
plot_circlize(circ_data,do.label = T, pt.size = 0.1,contour.levels = c(0.5, 0.8),col.use = col,
              bg.color = 'white', kde2d.n = 1000, repel = T, label.cex = 1.5)
add_track(circ_data, group = "stage",colors = RColorBrewer::brewer.pal(9, "Set1"), 
          track_num = 2) ## can change it to one of the columns in the meta data of your seurat object
add_track(circ_data, group = "author", colors = c(brewer.pal(12, "Set3"),   # 高区分度，适合分类
                                                  brewer.pal(8, "Dark2"),   # 深色系，色盲友好
                                                  brewer.pal(9, "Pastel1")[1:6]), 
          track_num = 3) ## can change it to one of the columns in the meta data of your seurat object
add_track(circ_data, group = "donor_id",colors = dittoColors(), 
          track_num = 4)
add_track(circ_data, group = "tissue",colors = tissue_colors, 
          track_num = 5)

dev.off()
