
knitr::opts_chunk$set(echo = TRUE)
library(Seurat)
library(SeuratObject)
library(Matrix)
library(irlba)
library(ggsci)
mycol <- c(pal_d3()(7),pal_aaas()(7),pal_uchicago()(7),pal_jama()(7))

work_dir <- "D:/ProgramData/Glioma-immunelncRNA/TAMs-MES"
setwd(work_dir)
getwd()
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

color_list <- ggplotColours(n=6)
color_list <- c(color_list, "grey")
gbm.combined_subset<-readRDS('D:/ProgramData/Spatial Gene Expression/Macrophage-mediated myelin recycling fuels brain cancer malignancy/GEO Data processed/scrna-seq/gbm-complete-subset.rds')
DefaultAssay(gbm.combined_subset) <- "RNA"

# Manually annotated cell types based on cell specific markers


## Assay is set to "RNA" as we will compare RNA expression between clusters
DefaultAssay(gbm.combined_subset) <- "RNA"
p3 <- DotPlot(gbm.combined_subset, group.by = "seurat_clusters", features = c("P2ry12", "Cd74", "Olig1" , "Gfap", "Pecam1", "Cd3e", "Itga4", "Mki67", "Des", "Ccr7", "Clec9a"))

p3

DimPlot(gbm.combined_subset, reduction = "umap", label = F)
# Split dataset for upcoming analyses and figures

## Dataset is split into a macrophage and tumor cell dataset
gbm.myeloid <- gbm.combined_subset[, gbm.combined_subset$celltype %in% c("MG", "MDMs")]

DefaultAssay(gbm.myeloid) <- "integrated"
gbm.myeloid <- ScaleData(gbm.myeloid, verbose = FALSE)
gbm.myeloid <- RunPCA(gbm.myeloid, npcs = 30, verbose = FALSE)

gbm.myeloid <- RunUMAP(gbm.myeloid, reduction = "pca", dims = 1:25)
gbm.myeloid <- FindNeighbors(gbm.myeloid, reduction = "pca", dims = 1:25)
gbm.myeloid <- FindClusters(gbm.myeloid, resolution = 0.26, algorithm = 1)

## Prepping tumor cells
gbm.tumor  <- gbm.combined_subset[, gbm.combined_subset$celltype %in% c("Tumor cells")]

DefaultAssay(gbm.tumor) <- "integrated"
gbm.tumor <- ScaleData(gbm.tumor, verbose = FALSE)
gbm.tumor <- RunPCA(gbm.tumor, npcs = 30, verbose = FALSE)

gbm.tumor <- RunUMAP(gbm.tumor, reduction = "pca", dims = 1:30)
gbm.tumor <- FindNeighbors(gbm.tumor, reduction = "pca", dims = 1:30)
gbm.tumor <- FindClusters(gbm.tumor, resolution = 0.3, algorithm = 1)
# Saving data for downstream analyses

save(gbm.myeloid,file="gbm-complete-myeloid.Rda")
save(gbm.tumor,file="gbm-complete-tumor.Rda")



knitr::opts_chunk$set(echo = TRUE)

# Figure 2B visium
library(Seurat)
library(nichenetr)
library(readxl)
library(dplyr)
library(readr)
library(nichenetr)
library(dplyr)
library(hdf5r)
library(ggplot2)


# Set working directory to folder "pre-processed visium" 
gbm.merged<-readRDS('D:/ProgramData/Spatial Gene Expression/Macrophage-mediated myelin recycling fuels brain cancer malignancy/GEO Data processed/visium/visium_merged_complete.rds')

table(gbm.merged$tissue)
gbm.merged@images
SpatialDimPlot(gbm.merged,ncol = 3)


SpatialDimPlot(gbm.merged, images = c( "Ink4a_Prim_S3"), group.by = "Location" ,  stroke = 0,   image.alpha = 0, alpha = 1) 
SpatialDimPlot(gbm.merged, images = c("Ink4a_Prim_S3"), group.by = "Cell.Subtype",  stroke = 0,   image.alpha = 0, alpha = 1) 

SpatialDimPlot(gbm.merged,  stroke = 0,   image.alpha = 0, alpha = 1, images = c("Ink4a_Prim_S3"), cells.highlight = CellsByIdentities(gbm.merged), facet.highlight = TRUE,  ncol = 4)
SpatialDimPlot(gbm.merged,  stroke = 0,   image.alpha = 0, alpha = 1, images = c("Ink4a_Prim_S3"), cells.highlight = CellsByIdentities(gbm.merged), facet.highlight = TRUE,  ncol = 4)


SpatialDimPlot(gbm.merged, group.by = "Location", images = c("Ink4a_Prim_S3")) + scale_fill_manual(values = c("orange", "violetred1", "darkolivegreen3","lightblue", "grey")) 
SpatialDimPlot(gbm.merged, group.by = "Cell.Subtype", images = c("Ink4a_Prim_S3"))


SpatialDimPlot(gbm.merged, group.by = "TME",   stroke = 0,   image.alpha = 0, alpha = 1, 
               images = c("Ink4a_Prim_S3")) + scale_fill_manual(values = color_list)

SpatialDimPlot(gbm.merged,  stroke = 0,   image.alpha = 0, alpha = 1, images = c("Ink4a_Prim_S3"), cells.highlight = CellsByIdentities(gbm.merged), facet.highlight = TRUE,  ncol = 4)

## Figure 2E-G
SpatialDimPlot(gbm.merged[, gbm.merged$TME %in% c("MG1-P2RY12",	"MG2-TNF", "MG3-GPNMB" ,	"MDM1-CCR7",	"MDM2-H2-EB1", "MDM3-GPNMB")], group.by = "GPNMBhigh_sig",   stroke = 0,   image.alpha = 1, alpha = 1, 
               images = c("Ink4a_Prim_S3")) + scale_fill_manual(values =  c("#FFA07A", "#AB82FF"))
SpatialDimPlot(gbm.merged, group.by = "Cell.Subtype",stroke = 0,  images = c("Ink4a_Prim_S3"), image.alpha = 0) + scale_fill_manual(values = c("darkolivegreen3", "pink", "tomato3","darkgoldenrod1", "darkgoldenrod3", "royalblue3", "NA" ))
SpatialDimPlot(gbm.merged, group.by = "Location",stroke = 0,  images = c("Ink4a_Prim_S3"), image.alpha =  0) + scale_fill_manual(values = c("orange", "violetred1", "darkolivegreen3","lightblue", "grey")) 

## Supplementary Figure 7C-E (additional representitive images of primary tumor)
SpatialDimPlot(gbm.merged[, gbm.merged$TME %in% c("MG1-P2RY12",	"MG2-TNF", "MG3-GPNMB" ,	"MDM1-CCR7",	"MDM2-H2-EB1", "MDM3-GPNMB")], group.by = "GPNMBhigh_sig",   stroke = 0,   image.alpha = 1, alpha = 1, 
               images = c("Ink4a_Prim_S3")) + scale_fill_manual(values =  c("#FFA07A", "#AB82FF"))
SpatialDimPlot(gbm.merged, group.by = "Cell.Subtype",stroke = 0,  images = c("Ink4a_Prim_S8_2"), image.alpha = 0) + scale_fill_manual(values = c("darkolivegreen3", "pink", "tomato3","darkgoldenrod1", "darkgoldenrod3", "royalblue3", "NA" ))
SpatialDimPlot(gbm.merged, group.by = "Location",stroke = 0,  images = c("Ink4a_Prim_S8_2"), image.alpha =  0) + scale_fill_manual(values = c("orange", "violetred1", "darkolivegreen3","lightblue", "grey")) 


######gbm.stRNAgbm.stRNAgbm.stRNAgbm.stRNAgbm.stRNAgbm.stRNAgbm.stRNAgbm.stRNAgbm.stRNAgbm.stRNAgbm.stRNAgbm.stRNAgbm.stRNAgbm.stRNAgbm.stRNAgbm.stRNA
gbm.stRNA<-readRDS('D:/ProgramData/Glioma-immunelncRNA/spatial-intergrate-gbmap/stRNA-anno-level234.rds')
table(gbm.stRNA$cna_bin)
gbm.stRNA$Metaprograms<-gbm.stRNA$cna_bin
library(forcats)
gbm.stRNA$Metaprograms<- fct_recode(gbm.stRNA$Metaprograms,
                                "Inflammatory-Mac" = "non_malignant",
                                "Prolif-Metab" = "mix_low",
                                "MES-hypoxia" = "mix_high",
                                "Phagocytic suppressive" = "malignant"
)
table(gbm.stRNA$Metaprograms)
Idents(gbm.stRNA) <- "Metaprograms"

gbm.stRNA@images
SpatialDimPlot(gbm.stRNA,ncol = 5)
colnames(gbm.stRNA@meta.data)


SpatialDimPlot(gbm.stRNA, images = c("UKF269"), group.by = "Metaprograms",  stroke = 0,   image.alpha = 0, alpha = 1) 

SpatialDimPlot(gbm.stRNA,  stroke = 0,   image.alpha = 0, alpha = 1, images = c("UKF269"), cells.highlight = CellsByIdentities(gbm.stRNA), facet.highlight = TRUE,  ncol = 4)


SpatialDimPlot(gbm.stRNA, group.by = "ivygap", images = c("UKF269")) + scale_fill_manual(values = c("orange", "violetred1", "darkolivegreen3","lightblue", "grey")) 
SpatialDimPlot(gbm.stRNA, group.by = "L3_second_type", images = c("UKF269"),image.alpha = 0)+ scale_fill_manual(values =  mycol)
SpatialDimPlot(gbm.stRNA, group.by = "L2_second_type", images = c("UKF269"))

SpatialDimPlot(gbm.stRNA, group.by = "L2_second_type",   stroke = 0,   image.alpha = 0, alpha = 1,images = c("UKF269")) + scale_fill_manual(values =  mycol)
SpatialDimPlot(gbm.stRNA, group.by = "L2_first_type",   stroke = 0,   image.alpha = 0, alpha = 1,images = c("UKF269")) 
SpatialDimPlot(gbm.stRNA, group.by = "layer",   stroke = 0,   image.alpha = 0, alpha = 1,images = c("UKF269")) + scale_fill_manual(values =  mycol)

SpatialDimPlot(gbm.stRNA,  stroke = 0,   image.alpha = 0, alpha = 1, images = c("UKF269"), cells.highlight = CellsByIdentities(gbm.stRNA), facet.highlight = TRUE,  ncol = 4)

## Figure 2E-G
SpatialDimPlot(gbm.stRNA[, gbm.stRNA$L2_second_type %in% c("Myeloid")], group.by = "L3_second_type",   stroke = 0,   image.alpha = 0, alpha = 1, 
               images = c("UKF269")) + scale_fill_manual(values =  mycol)

SpatialDimPlot(gbm.stRNA, group.by = "ivygap",stroke = 0,  images = c("UKF269"), image.alpha = 0) + scale_fill_manual(values = c("darkolivegreen3", "pink", "tomato3","darkgoldenrod1", "darkgoldenrod3", "royalblue3", "NA" ))
SpatialDimPlot(gbm.stRNA, group.by = "Metaprograms",stroke = 0,  images = c("UKF269"), image.alpha =  0) + scale_fill_manual(values = c("orange", "violetred1", "darkolivegreen3","lightblue", "grey")) 

## Supplementary Figure 7C-E (additional representitive images of primary tumor)
SpatialDimPlot(gbm.stRNA[, gbm.stRNA$TME %in% c("MG1-P2RY12",	"MG2-TNF", "MG3-GPNMB" ,	"MDM1-CCR7",	"MDM2-H2-EB1", "MDM3-GPNMB")], group.by = "GPNMBhigh_sig",   stroke = 0,   image.alpha = 1, alpha = 1, 
               images = c("UKF269")) + scale_fill_manual(values =  c("#FFA07A", "#AB82FF"))
SpatialDimPlot(gbm.stRNA, group.by = "Cell.Subtype",stroke = 0,  images = c("Ink4a_Prim_S8_2"), image.alpha = 0) + scale_fill_manual(values = c("darkolivegreen3", "pink", "tomato3","darkgoldenrod1", "darkgoldenrod3", "royalblue3", "NA" ))
SpatialDimPlot(gbm.stRNA, group.by = "Location",stroke = 0,  images = c("Ink4a_Prim_S8_2"), image.alpha =  0) + scale_fill_manual(values = c("orange", "violetred1", "darkolivegreen3","lightblue", "grey")) 


# 提取子集
gbm_UKF269 <- subset(gbm.stRNA, subset = orig.ident == "UKF269")
gbm_UKF269@images <- list(UKF269 = gbm_UKF269@images$UKF269)# 现在应该只有1个图像
length(gbm_UKF269@images)  
# 可视化（无背景H&E图）
SpatialDimPlot(
  gbm_UKF269,
  image.alpha = 0,        # 移除背景
  stroke = NA,             # 去除点描边
  pt.size.factor = 1.5     # 调整点大小
) + theme_void()
save(gbm_UKF269,file = 'gbmUKF269.rds')


DimPlot(gbm.stRNA, reduction = "umap", label = TRUE)|DimPlot(gbm.stRNA, reduction = "tsne", label = TRUE)



#可视化每个sport的空间计数
SpatialFeaturePlot(gbm_UKF269, features = "nCount_Spatial")


##SCT标准化
gbm_UKF269 <- SCTransform(gbm_UKF269, assay = "Spatial", verbose = FALSE)
gbm_UKF269 <- RunPCA(gbm_UKF269, assay = "SCT", verbose = FALSE) 

##数据聚类
gbm_UKF269<- FindNeighbors(gbm_UKF269, reduction = "pca", dims = 1:10)
gbm_UKF269 <- FindClusters(gbm_UKF269, verbose = FALSE,resolution = 0.6)
p1<-SpatialPlot(gbm_UKF269, label = TRUE, label.size = 5)

#UMAP降维
gbm_UKF269 <- RunUMAP(gbm_UKF269, reduction = "pca", dims = 1:10)
p2 <- DimPlot(gbm_UKF269, reduction = "umap", label = TRUE)
p1+p2

SpatialFeaturePlot(gbm_UKF269, features =c("SOX10","SOD2"))
SpatialPlot(gbm_UKF269,label=TRUE,group.by = 'seurat_clusters',label.size=8)
save(gbm_UKF269,file = 'gbm_UKF269S1.rdata')

# 
# 
# ##实战2：空间区域划分及差异分析
# ##加载R包
# library(Seurat)
# library(ggplot2)
# library(tidyverse)
# 
# ##载入实战1保存的Seurat对象
# load('gbm_UKF269S1.rdata')
# 
# #视化Cluster分布及H&E组织切片
# p1<-SpatialPlot(gbm_UKF269, label = TRUE, label.size = 5)
# p2<-SpatialPlot(gbm_UKF269,pt.size.factor = 0.6)+NoLegend()
# p1+p2
# 
# #区域定义
# gbm_UKF269@meta.data$Region<-NA
# gbm_UKF269@meta.data$Region[gbm_UKF269@meta.data$seurat_clusters %in% c('0','3','5')] <- "Normal"
# gbm_UKF269@meta.data$Region[gbm_UKF269@meta.data$seurat_clusters %in% c('4')] <- "Tumor"
# gbm_UKF269@meta.data$Region[gbm_UKF269@meta.data$seurat_clusters %in% c('1','2','7','6')] <- "Tumor core"
# SpatialPlot(gbm_UKF269, label = TRUE, label.size = 5,group.by = 'Region',cols = c('Normal'='#4b5cc4','Transition'='#FE8D3C','Tumor'='#AA0000'))

##区域间差异分析
#切换Idents为上面定义的Metaprograms
Idents(gbm_UKF269)<-gbm_UKF269$Metaprograms

#找到各区域的标记物，并只报道阳性位点
markers <- FindAllMarkers(gbm_UKF269, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#找到各区域top10基因
top10<-markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)


#可视化top10基因热图
DoHeatmap(gbm_UKF269, features = top10$gene) + NoLegend()


library(clusterProfiler)

#将基因SYMBOL转换为ENTREZID
gid <- bitr(unique(markers$gene), 'SYMBOL', 'ENTREZID', OrgDb= 'org.Hs.eg.db')
markers <- full_join(markers, gid, by=c('gene' = 'SYMBOL'))

#KEGG通路富集分析
KEGG = compareCluster(ENTREZID ~ cluster, data = markers, fun='enrichKEGG')

#可视化各区域TOP5通路结果
dotplot(KEGG, label_format=40) + theme(axis.text.x = element_text(angle=45, hjust=1)) + scale_color_gradient(high="#4b5cc4",low="#FE8D3C")

#保存分析结果，以便下期使用
save(gbm_UKF269,file = 'gbm_UKF269S2.rdata')






###实战6:空间分化轨迹推断
###基于Monocle2
library(Seurat)
library(monocle)
library(tidyverse)
library(patchwork)

#载入实战2保存的 gbm_UKF269 空转数据
load('gbm_UKF269S2.rdata')
Idents(gbm_UKF269)<-"Metaprograms" 
#可视化空间区域
SpatialPlot(gbm_UKF269,label=TRUE,
            label.size=5)

##提取出Transition与Tumor的spot进行后续分析
subGBM<-gbm_UKF269
#可视化区域
SpatialPlot(subGBM,label=TRUE,
            label.size=5,
            group.by = "Metaprograms")

#可视化clusters
SpatialPlot(subGBM,label=TRUE,
            group.by = 'seurat_clusters',
            label.size=8)
##获取原始表达矩阵
data<-as(as.matrix(subGBM@assays$Spatial@counts),'sparseMatrix')

#构建mycds对象
pd<-new('AnnotatedDataFrame',data =subGBM@meta.data)
fData<-data.frame(gene_short_name =row.names(data),row.names=row.names(data))
fd<-new('AnnotatedDataFrame',data =fData )

mycds<-newCellDataSet(data,
                      phenoData=pd,
                      featureData=fd,
                      expressionFamily=negbinomial.size())#从seurat到monocle

##估计size factor和离散度
mycds<-estimateSizeFactors(mycds)
mycds<-estimateDispersions(mycds,cores=1,relative_expr = TRUE)

##使用monocle选择高变基因
disp_table<-dispersionTable(mycds)
#过滤低质量的细胞
disp.genes<-subset(disp_table,mean_expression>= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
##轨迹构建基因可视化
mycds <- setOrderingFilter(mycds,disp.genes)
plot_ordering_genes(mycds)
##降维
mycds <- reduceDimension(mycds, max_components = 2,method ='DDRTree')
##默认排序
##服务器来吧
saveRDS(mycds,file = 'reduceDimensionmycds.RDS')
mycds<-orderCells(mycds)##服务器了
mycds<-readRDS('orderCellscds-UKF269.RDS')
##Pseudotime轨迹可视化
p1<-plot_cell_trajectory(mycds, color_by="Pseudotime")
p2<-plot_cell_trajectory(mycds, color_by="seurat_clusters")
p3<-plot_cell_trajectory(mycds, color_by="State")
p1+p2+p3

#把细胞的分化时间映射到组织切片上
##提取伪时间数据
spot_Pseudotime <- data.frame(pData(mycds)$Pseudotime)
rownames(spot_Pseudotime) <- rownames(subGBM@meta.data)

##伪时间数据加至gbm_UKF269@meta.data
gbm_UKF269[['Pseudotime']] <- 0
gbm_UKF269[['Pseudotime']][rownames(spot_Pseudotime),] <- spot_Pseudotime

##伪时间空间可视化
SpatialFeaturePlot(gbm_UKF269,features=c("Pseudotime"))

##寻找拟时差异基因
Time_diff <- differentialGeneTest(mycds[disp.genes,], cores = 8, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_genes <- Time_diff%>% pull(gene_short_name) %>% as.character()

#num_clusters=4,可根据热图效果自行调整
p<-plot_pseudotime_heatmap(mycds[Time_genes,], num_clusters=4, show_rownames=F, return_heatmap=T)



##获取TOP5显著差异基因
top_genes <- Time_diff %>%
  slice_max(order_by = -qval, n = 5) %>%
  pull(gene_short_name) %>%
  as.character()
top_genes

SpatialFeaturePlot(subGBM,features=c("PLA2G2A","S100A10","CHI3L1","TMSB10","TF", "SPP1" ))

SpatialFeaturePlot(subGBM,features=c("GOT1","GLS2" ,"GOT2"
                                     ,"GPX4","HMOX1","TMEM164","AIFM2",
                                     "NFE2L2","SLC39A7","SLC7A11","ITGAM"))

##基因随细胞状态轨迹图
cds_subset <- mycds[top_genes[1:6],]
p1<-plot_genes_in_pseudotime(cds_subset,color_by = "Pseudotime")+
  scale_color_viridis_c(option = "viridis")
p2<-plot_genes_in_pseudotime(cds_subset,color_by = "State")
p1+p2


#上述热图设置为4个cluster,这里的K值则保持一致
clusters <- cutree(p$tree_row, k = 4)
geneCluster <- data.frame(clusters)
geneCluster[,1] <- as.character(geneCluster[,1])
colnames(geneCluster) <- "Clusters"
geneCluster$gene=rownames(geneCluster)
table(geneCluster$Clusters)

## Cluster基因_GO富集分析
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(RColorBrewer)

#转换id
gid <- bitr(unique(geneCluster$gene), 'SYMBOL', 'ENTREZID', OrgDb= 'org.Hs.eg.db')  # 人 Hs 鼠 Mm
Genelist <- full_join(geneCluster, gid, by=c('gene' = 'SYMBOL'))

##GO:BP 富集分析
GO_BP = compareCluster(ENTREZID ~ Clusters, 
                       data = Genelist, 
                       fun='enrichGO',
                       OrgDb = 'org.Hs.eg.db',
                       ont = "BP") 
dotplot(GO_BP,label_format = 60)





























