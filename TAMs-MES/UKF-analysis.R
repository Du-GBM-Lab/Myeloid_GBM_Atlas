###实战6:空间分化轨迹推断
###基于Monocle2
library(Seurat)
library(monocle)
library(tidyverse)
library(patchwork)
library(qs)

gbm_UKF269<-qread('gbmUKF269.qs')
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
p2<-plot_cell_trajectory(mycds, color_by="Metaprograms")
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
SpatialFeaturePlot(subGBM,features=c("CAMK2N1","C1QB","ID3","HMGN2","FABP3", "SPP1" ))

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



##Seurat转换为SPATA2对象
#Convert
library(SPATA2)
names(subGBM@assays)

# use Spatial assay 
spata_object1 <-
  asSPATA2(
    object = subGBM,
    sample_name = "subGBM",
    platform = "VisiumSmall",
    img_scale_fct = "lowres",
    assay_name = "Spatial", 
    assay_modality = "gene"
  )

show(spata_object1)

# seurat_clusters
plotSurface(object = spata_object1, color_by = "seurat_clusters")
##运行 CNV 分析
spata_object1 <-
  runCNV(
    object = spata_object1,
    #设置结果输出目
    directory_cnv_folder = "./UKF269-SPATA/", 
    cnv_prefix = "Chr"
  )
#CNV热图
plotCnvHeatmap(object = spata_object1, across = "Metaprograms", clrp = "npg")
qsave(spata_object1,file = 'spata_object1.qs')
spata_object1<-qread('spata_object1.qs')
##从热图结果看，Cluster1和4的CNV更高
#指定亚群，SNAP和chrom可视化
plotCnvHeatmap(
  object = spata_object1, 
  across = "seurat_clusters", 
  across_subset = c('1','4'), # dont show the transition part
  meta_vars = "SNAP29",  # visualize SNAP29 expression on the left
  meta_vars_clrs = c("SNAP29" = "inferno"), # with the inferno color spectrum 
  chrom_subset = c("6", "7", "11", "12", "19"), # only show these chromosomes
  ggpLayers = list(arm = list(legendNone())) # remove the chrom arm legend
)
plotSurface(spata_object1, color_by = "SNAP29", pt_clrsp = "inferno")
##折线图可视化
plotCnvLineplot(
  object = spata_object1,
  across = "Metaprograms", 
  n_bins_bcsp = 1000,
  n_bins_genes = 1000,
  nrow = 4
)
##空间可视化染色体改变
#可视化空间染色体改变 
getCnvFeatureNames(object = spata_object1) %>% head()

# are part of all feature names
getFeatureNames(object = spata_object1) %>% head()

plotSurface(
  object = spata_object1, 
  color_by = "Chr7", 
  pt_clrsp = "Reds"
)

plotSurface(
  object = spata_object1, 
  color_by = "Chr10", 
  pt_clrsp = "Oslo" 
)

save(spata_object1,file = 'U==KF269_SPATA2.rdata')

###空间共定位分析-MISTy
library(mistyR)
library(CARD)
library(Seurat)
library(mistyR)
library(distances)
library(future)
library(tidyverse) 
library(recipes)

gbm_UKF269<-qread('gbmUKF269.qs')
SpatialPlot(gbm_UKF269)



##共表达网络分析-hdWGCNA

library(Seurat)
library(hdWGCNA)
library(WGCNA)
library(tidyverse)
library(cowplot)
library(patchwork)

#载入实战2保存的结果

UKF269T<-qread('gbmUKF269.qs')
SpatialPlot(UKF269T)
colnames(UKF269T@meta.data)


##启用网络分析的并行处理（可选）
enableWGCNAThreads(nThreads = 16)

#参考教程，重新走标准化流程，没有用SCT方法
UKF269T <- UKF269T %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

#clustering and umap
UKF269T <- FindNeighbors(UKF269T, dims = 1:30)
UKF269T <- FindClusters(UKF269T,verbose = TRUE)
UKF269T <- RunUMAP(UKF269T, dims = 1:30)

DimPlot(UKF269T, label=TRUE, reduction = "umap", group.by = "seurat_clusters") 
SpatialDimPlot(UKF269T, label = TRUE, label.size = 3, group.by = "Metaprograms")
#提取位置信息，并加入到UKF269T@meta.data
coord<-UKF269T@images[["UKF269"]]@coordinates
UKF269T<-AddMetaData(UKF269T,metadata = coord)
UKF269T@meta.data[1:5,1:12]
#构建 metaspots
UKF269T <- SetupForWGCNA(
  UKF269T,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = "vis"
)


UKF269T <- MetaspotsByGroups(
  UKF269T,
  group.by = c("Metaprograms"),
  ident.group = "Metaprograms",
  assay = 'Spatial'
)

UKF269T  <- NormalizeMetacells(UKF269T)


#设置表达式矩阵，将group.by和group_name设置为NULL以包含所有点
UKF269T  <- SetDatExpr(
  UKF269T,
  assay='Spatial',
  group.by=NULL,
  group_name = NULL
)

#测试不同的软阈值
UKF269T <- TestSoftPowers(UKF269T)
plot_list <- PlotSoftPowers(UKF269T)

wrap_plots(plot_list, ncol=2)

#构建共表达网络
UKF269T <- ConstructNetwork(
  UKF269T,
  tom_name='test',
  overwrite_tom=TRUE
)

#绘制树状图
PlotDendrogram(UKF269T, main='Spatial hdWGCNA dendrogram')

UKF269T <- ModuleEigengenes(UKF269T)
UKF269T <- ModuleConnectivity(UKF269T)

#重置带有前缀 “SM” 的模块名称(可选)
UKF269T <- ResetModuleNames(
  UKF269T,
  new_name = "M"
)

#获取模块特征基因和基因模块分配表
MEs <- GetMEs(UKF269T)
modules <- GetModules(UKF269T)

#去除grey模块
mods <- levels(modules$module); mods <- mods[mods != 'grey']

#将ME添加到seurat元数据中，这样我们就可以用seurat函数绘制它
UKF269T@meta.data <- cbind(UKF269T@meta.data, MEs)

#使用Seurat的DotPlot函数绘图
p <- DotPlot(UKF269T, features=mods, group.by = 'Metaprograms', dot.min=0.1)

#翻转x/y轴，旋转轴标签，并更改配色方案：
p <- p +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  xlab('') + ylab('')

p
#直接在空间坐标上可视化模块
SpatialFeaturePlot(
  UKF269T,
  features = mods,
  alpha = c(0.1, 1),
  ncol = 3
)


#在共表达网络上执行UMAP嵌入
UKF269T <- RunModuleUMAP(
  UKF269T,
  n_hubs = 5,
  n_neighbors=15,
  min_dist=0.3,
  spread=1
)

#绘制网络图
ModuleUMAPPlot(
  UKF269T,
  edge.alpha=0.5,
  sample_edges=TRUE,
  keep_grey_edges=FALSE,
  edge_prop=0.075, 
  label_hubs=5 
)
library(igraph)

#默认输出各模块TOP25个hub基因，图片保存在指定的输出路径下
ModuleNetworkPlot(
  UKF269T,
  outdir = './ModuleNetworks/'
)
#保存结果
save(UKF269T,file = 'UKF269T_hdWGCNA.rdata')

#查看各模块基因数量
table(modules$module)
