# Prepare Data In the first step we will load the data and prepare the multi-omic methylation and RNA-seq information:
setwd("/home/data/t060202/Spatial-GSVA-duijian/epigenetic in glioma/")
setwd("/home/data/t060202/Spatial-GSVA-duijian/Spatial_Glioma-CELL/intergreted/")

.libPaths(c('~/R/x86_64-pc-linux-gnu-library/4.3/'))


library(harmony)
library(WGCNA)
library(hdWGCNA)
# enable parallel processing for network analysis (optional)
enableWGCNAThreads(nThreads = 20)
# using the cowplot theme for ggplot
library(ggplot2)
set.seed(12345)
library(SPATA2)
library(tidyverse)
library(Seurat) ####4
library(ggsci)
library(ggplot2)
library(bigmemory)
library(doMC)
library(patchwork)
library(future)
library(stringr)
library(dplyr)
library(HGNChelper)
library(future)
theme_set(theme_cowplot())
mycol <- c(pal_d3()(7),pal_aaas()(7),pal_uchicago()(7),pal_jama()(7))
plan("multisession", workers = 24)
plan()
#设置可用的内存
options(future.globals.maxSize = 100000* 1024^3)
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

p_mt.list <- lapply(spatial.list, function(x){
  p1 <- VlnPlot(x, features = "percent.mt") + ggtitle(unique(x$orig.ident)) +
    theme(legend.position = "none", axis.text.x = element_blank())
  p2 <- SpatialFeaturePlot(x, features = "percent.mt") + theme(legend.position = "right")
  p <- p1|p2
  p })
wrap_plots(p_mt.list, ncol = 5)
# 整合样本
stRNA <- merge(spatial.list[[1]], spatial.list[2:length(spatial.list)])
# make a dataframe containing the image coordinates for each sample
image_df <- do.call(rbind, lapply(names(stRNA@images), function(x){
  stRNA@images[[x]]@coordinates
}))

# merge the image_df with the Seurat metadata
new_meta <- merge(stRNA@meta.data, image_df, by='row.names')

# fix the row ordering to match the original seurat seurat_vhd
rownames(new_meta) <- new_meta$Row.names
ix <- match(as.character(colnames(stRNA)), as.character(rownames(new_meta)))
new_meta <- new_meta[ix,]

# add the new metadata to the seurat seurat_vhd
stRNA@meta.data <- new_meta

head(image_df)
#stRNA = JoinLayers(stRNA)
#stRNA <-stRNA[,unname(which(colSums(GetAssayData(stRNA,,assay="Spatial"))!=0))]
# 标准化
# normalization, feature selection, scaling, and PCA
stRNA <- stRNA %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

# Louvain clustering and umap
stRNA <- FindNeighbors(stRNA, dims = 1:30)
stRNA <- FindClusters(stRNA,verbose = TRUE)
stRNA <- RunUMAP(stRNA, dims = 1:30)


# show the UMAP
p1 <- DimPlot(stRNA, label=TRUE, reduction = "umap", group.by = "seurat_clusters") + NoLegend()
p1
p2 <- SpatialDimPlot(stRNA,label = T, alpha = c(0.3, 1), 
                     facet.highlight = TRUE, ncol = 5,
                     pt.size.factor = 2,label.size = 3.5)
p2
common_cells <- intersect(colnames(stRNA), rownames(metadata))
dim(stRNA)
stRNA_filtered <- subset(stRNA, cells = common_cells)
dim(stRNA_filtered)
dim(metadata)
seurat_obj <- AddMetaData(
  object = stRNA_filtered,
  metadata = metadata
)
rm(stRNA_filtered)
head(seurat_obj@meta.data)
table(seurat_obj$sample)
print(identical(seurat_obj$orig.ident, seurat_obj$sample))  # 返回 TRUE/FALSE
DimPlot(seurat_obj, reduction = "umap", label = FALSE,group.by = 'sample',pt.size = 0.15,cols = mycol,label.size = 7)

saveRDS(seurat_obj, file="seurat_obj_before_StepWGCNA.Rds")
seurat_obj<-readRDS('stRNA-anno-level234.rds')

seurat_obj$ivygap <- factor(as.character(seurat_obj$ivygap), levels=c('CT', 'LE','MVP','PAN'))

###Construct metaspots
seurat_obj <- SetupForWGCNA(seurat_obj,gene_select = "fraction",fraction = 0.05,wgcna_name = "vis")
seurat_obj <- MetaspotsByGroups(seurat_obj,
                                group.by = c("orig.ident","ivygap"),
                                ident.group = "ivygap",assay = 'Spatial')
seurat_obj  <- NormalizeMetacells(seurat_obj)
m_obj <- GetMetacellseurat_vhd(seurat_obj)
m_obj
# set up the expression matrix, set group.by and group_name to NULL to include all spots
seurat_obj  <- SetDatExpr(seurat_obj,group.by=NULL,group_name = NULL)
# test different soft power thresholds
seurat_obj <- TestSoftPowers(seurat_obj)
wrap_plots(PlotSoftPowers(seurat_obj), ncol=2)
saveRDS(seurat_obj,file = 'hdWGCNA_TestSoftPowers_seurat.obj-st.rds')
seurat_obj<-readRDS('hdWGCNA_TestSoftPowers_seurat.obj-st.rds')
# construct co-expression network:
seurat_obj <- ConstructNetwork(
  seurat_obj,
  tom_name='test',
  overwrite_tom=TRUE
)

# plot the dendrogram
PlotDendrogram(seurat_obj, main='Spatial hdWGCNA dendrogram')
seurat_obj <- ModuleEigengenes(seurat_obj)
seurat_obj <- ModuleConnectivity(seurat_obj)
###ssgsea_seurat_obj########
library(irGSEA)
setwd("/home/data/t060202/Spatial-GSVA-duijian/Spatial_Glioma-CELL/intergreted/")
seurat_obj_irssgsea<- irGSEA.score(object = seurat_obj, assay = "Spatial", slot = "data", seeds = 123, ncores = 24, 
                                   min.cells = 3, min.feature = 0, custom = F, 
                                   geneset = NULL, msigdb = T, species = "Homo sapiens", category = "H", 
                                   subcategory = NULL, geneid = "symbol", method = c( "ssgsea"), 
                                   aucell.MaxRank = NULL, ucell.MaxRank = NULL, kcdf = 'Gaussian')
ssgsea_scores <- seurat_obj_irssgsea@assays$ssgsea@scale.data
gene_sets <- rownames(ssgsea_scores)
for (gs in gene_sets) {
  seurat_obj <- AddMetaData(seurat_obj, metadata = ssgsea_scores[gs, ], col.name = gs)
}
seurat.obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 5,
  n_neighbors=5,
  min_dist=0.1,
  spread=2,
  wgcna_name = 'vis',
  target_weight=0.05,
  supervised=TRUE
)

UMAP=ModuleUMAPPlot(
  seurat.obj,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.2, 
  label_hubs=2,
  return_graph=T,
  keep_grey_edges=FALSE)
### Save Data

saveRDS(list(objet=seurat.obj, UMAP=UMAP), file="Analysis_StepWGCNA.Rds")

dlist<-readRDS("Analysis_StepWGCNA.Rds")
# 提取 Seurat 对象
seurat.obj <- dlist$objet  # 或 data[["objet"]] 
# 提取 UMAP 结果
umap_result <- dlist$UMAP  # 或 data[["UMAP"]]
## Plot Gene Expression UMAP
library(scattermore)
library(igraph)
table(seurat.obj@misc$vis$module_umap$color)
plotdf <- seurat.obj@misc$vis$module_umap %>% filter(color!="grey60")
plotdf2 <- as_long_data_frame(UMAP)
plotdf2 <- plotdf2 %>% filter(from_name %in% plotdf$gene)

# Create Plot of Expression Network
plot_p <- ggplot()+theme_void()+coord_fixed()
dat1 <- plotdf2 %>% filter(from_kME>0.5)
plot_p <- 
  plot_p+
  geom_segment(data=dat1, 
               mapping = aes(x=from_UMAP1, y=from_UMAP2, 
                             xend=to_UMAP1, yend=to_UMAP2), 
               size=dat1$from_kME*0.1, alpha=dat1$from_kME*0.1)

plot_p <- plot_p+geom_scattermore(data=plotdf, 
                                  mapping=aes(x=UMAP1, y=UMAP2),
                                  size=plotdf$kME*2, alpha=0.1, 
                                  color=gplots::col2hex(plotdf$color), pointsize=3)

library(ggrepel)

plot_p <- plot_p+
  geom_text_repel(data=plotdf %>% group_by(module) %>% 
                    summarise(UMAP1=mean(UMAP1), UMAP2=mean(UMAP2)), mapping=aes(x=UMAP1, y=UMAP2,label=module))

#plot_p <- 
#  plot_p+
#  geom_text_repel(data=plotdf %>% group_by(module) %>% 
#                    top_n(., 2, wt=kME), 
#                  mapping=aes(x=UMAP1, y=UMAP2,label=gene), size=2)
plot_p+
  ggtitle("Dimensional reduction of modules")+
  xlab("UMAP 1")+ylab("UMAP 2")+
  theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", size=0.5),
          axis.text.x = element_text(colour="black",
                                     angle = 90, vjust = 0.5, hjust=1), 
          axis.text.y = element_text(colour="black"))

######### Differential expressed modules across neural high/low
DMEs_all <- FindAllDMEs(
  seurat.obj,
  group.by = 'ivygap',
  wgcna_name = 'vis'
)

head(DMEs_all)

sig.modules <- 
  DMEs_all %>% 
  filter(avg_log2FC>0.5 & !is.infinite(avg_log2FC)) %>% 
  group_by(group) %>% 
  top_n(., 50, wt=-log(p_val)) %>% 
  arrange(desc(avg_log2FC), .by_group = T) %>% 
  ungroup() %>% 
  mutate(module=as.factor(module))
col <- gplots::col2hex(sig.modules$module)
names(col) <- sig.modules$module
sig.modules$module <- factor(sig.modules$module, levels = sig.modules$module)

#barplot
ggplot(sig.modules, aes(fill=module, y=avg_log2FC, x=group)) + 
  geom_bar(stat="identity",position=position_dodge())+
  theme_classic()+
  scale_fill_manual(values=col)+
  ggtitle("Differentially Expression")+
  xlab("")+ylab("")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1),
        axis.text.x = element_text(colour="black",
                                   angle = 90, vjust = 0.5, hjust=1), 
        axis.text.y = element_text(colour="black"))













plot_mod <- c("turquoise")
library(scattermore)

plotdf <- seurat.obj@misc$vis$module_umap %>% filter(color!="grey60") %>% filter(color %in% plot_mod)
plotdf2 <- as_long_data_frame(UMAP)
plotdf2 <- plotdf2 %>% filter(from_name %in% plotdf$gene)

# Create Plot of Expression Network
plot_p <- ggplot()+theme_void()+coord_fixed()
dat1 <- plotdf2 %>% filter(from_kME>0.5)


plot_p <- 
  plot_p+
  geom_segment(data=dat1, 
               mapping = aes(x=from_UMAP1, y=from_UMAP2, 
                             xend=to_UMAP1, yend=to_UMAP2), 
               size=dat1$from_kME*0.1, alpha=dat1$from_kME*0.3)

plot_p <- plot_p+geom_point(data=plotdf, 
                                  mapping=aes(x=UMAP1, y=UMAP2),
                                  size=plotdf$kME*2, alpha=0.1, 
                                  color=gplots::col2hex(plotdf$color))

plot_p <- plot_p+
  geom_text_repel(data=
                    plotdf %>% 
                    group_by(module) %>% 
                    summarise(UMAP1=mean(UMAP1), UMAP2=mean(UMAP2)), 
                  mapping=aes(x=UMAP1, y=UMAP2, label=module))

plot_p <- 
  plot_p+
  geom_text_repel(data=plotdf %>% group_by(module) %>% 
                    top_n(., 20, wt=kME), 
                  mapping=aes(x=UMAP1, y=UMAP2,label=gene),size=2, force=10, force_pull=10, max.overlaps=20,
                  segment.size=0.1)
plot_p+
  ggtitle("Dimensional reduction of modules")+
  xlab("UMAP 1")+ylab("UMAP 2")+
  theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", size=0.5),
          axis.text.x = element_text(colour="black",
                                     angle = 90, vjust = 0.5, hjust=1), 
          axis.text.y = element_text(colour="black"))

## Performe trait correlation
# Traits
seurat.obj@meta.data$ivygap <- as.factor(seurat.obj@meta.data$ivygap)
colnames(seurat.obj@meta.data)
traits <- c("ivygap","HALLMARK.HYPOXIA")

# list of traits to correlate
seurat.obj <- ModuleTraitCorrelation(seurat.obj,traits = traits)

seurat.obj <- ModuleTraitCorrelation(
  seurat.obj,
  traits = traits)

mt_cor <- GetModuleTraitCorrelation(seurat.obj)

## Remove modules without significants
keep <- 
  map(1:nrow(t(mt_cor$pval$all_cells < 0.05)), 
    ~any(t(mt_cor$pval$all_cells < 0.05)[.x,]==T)) %>% 
  unlist() %>% 
  which()

col=colorRampPalette(c(RColorBrewer::brewer.pal(9, "BrBG")))
library(ggcorrplot)
ggcorrplot(t(mt_cor$cor$all_cells)[keep, ] %>% as.data.frame(), 
           method = "circle", 
           outline.color = "black", 
           lab_size=1,
           insig="pch",
           pch="X",
           p.mat=t(mt_cor$pval$all_cells)[keep, ])+
  scale_fill_gradientn(colours = col(50), name="")+
  guides(fill = guide_colourbar(barwidth = 0.3, barheight = 8, ticks =F, frame.colour="black"), label=F)+
  ggtitle("Module correlation to treatment condition")+
  xlab("")+ylab("")+
  theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", size=1),
          axis.text.x = element_text(colour="black",
                                     angle = 90, vjust = 0.5, hjust=1), 
          axis.text.y = element_text(colour="black"))


## Gene Ontology Analysis:

modules <- c("blue", "yellow", "turquoise", "brown")

list_genes <- map(.x=modules, .f=function(i){
  genes <- 
    seurat.obj@misc$vis$wgcna_modules %>% 
    filter(module==i) %>% 
    rownames()
  
  genes <- c(unique(genes))
})

all <- unlist(list_genes)
names(list_genes) <- modules



library(DOSE)
library(enrichplot)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db) 
# 1. 去重处理
list_genes <- lapply(list_genes, function(genes) {
  mapIds(org.Hs.eg.db, keys = genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
})
list_genes <- lapply(list_genes, function(x) x[!is.na(x)])

# Run the Gene Ontology Analysis of the different clusters
CC <- compareCluster(list_genes, 
                     fun="enrichGO",
                     ont ="BP", 
                     keyType = "ENTREZID", 
                     OrgDb = "org.Hs.eg.db")
edo <- pairwise_termsim(CC)

## Plot the comparison
col <- colorRampPalette(c(RColorBrewer::brewer.pal(9, "Greens"))[4:9])

enrichplot::dotplot(edo, showCategory=3)+
  scale_color_gradientn(colours = col(50), name="")+
  guides(color = guide_colourbar(barwidth = 0.3, barheight = 8, 
                                 ticks =F, frame.colour="black"), label=F, size="none")+
  ggtitle("Gene Ontology")+
  xlab("")+ylab("")+
  theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", size=0.5),
          axis.text.x = element_text(colour="black",
                                     angle = 90, vjust = 0.5, hjust=0.5 ,size=5), 
          axis.text.y = element_text(colour="black", size=5))



## Get a Network representation of the enriched pathways in all cluster

#Set the random start point to get allways the same plots
set.seed(200)
plot <- 
  emapplot(edo, showCategory=50)



res_map <- 
  edo@compareClusterResult %>% 
  filter(Description %in% plot$data$name) %>% 
  group_by(Description, Cluster) %>% 
  summarise(p = mean(p.adjust)) %>% 
  ungroup() 

anno <- data.frame(name=unique(res_map$Description), type=map(.x=unique(res_map$Description), .f=function(i){
  c <- res_map %>% filter(Description==i) %>% pull(Cluster)
  return(c[res_map %>% filter(Description==i) %>% pull(p) %>% which.min()] %>% as.character())
  }) %>% unlist())

library(ggrepel)
plot$data <- plot$data %>% left_join(., anno, by="name")

col.test <- gplots::col2hex(unique(plot$data$type))
names(col.test) <- unique(plot$data$type)

ggplot(data=plot$data)+
  plot$layers[[1]]+
  geom_point(mapping = aes(x=x, y=y, fill=type, size=size), colour="black",pch=21)+
  scale_fill_manual(values=col.test)+
  scale_size(range=c(1, 2))+
  theme_bw() +
  guides(size="none")+
  xlab("Dim 1")+ylab("Dim 2")+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", size=0.5),
          axis.text.x = element_text(colour="black", size=0.5), 
          axis.text.y = element_text(colour="black", size=0.5))+
    coord_fixed()+
  geom_text_repel(mapping= aes(x,y, label=name),
                  size=2, force=10, force_pull=10, max.overlaps=20,
                  segment.size=0.1)+
  ggtitle("GeneSet Enrichment Analysis")
(".../SPATA_plot_extensions.R")
source(".../Run_extensions.R")
setwd("/home/data/t060202/Spatial-GSVA-duijian/Spatial_Glioma-CELL/intergreted/")

## New dataset with neurons
scRef <- readRDS("GBMAP_ssgsea.RDS")
colors <- readRDS("GBM_Neuron_colors.RDS")
table(colors$annotation_level_4)
table(scRef$annotation_level_2)
table(scRef$annotation_level_4)
colors$annotation_level_4[1:54] <- 
  str_replace_all(colors$annotation_level_4[1:54], "[ ]", "_") %>%
  str_replace_all(., "-", "_") %>%
  str_replace_all(., "/", "_")
col_cells <- colors$colors
names(col_cells) <- colors$annotation_level_4
## Get Plots of the modules:

# Plot the single cell data
f <- DimPlot(scRef, group.by  = "annotation_level_4")

library(scattermore)
f$data$annotation_level_4 <- factor(f$data$annotation_level_4, levels = colors$annotation_level_4)
p <- ggplot(data=f$data)+
  geom_scattermore(mapping=aes(UMAP_1,UMAP_2), pointsize = 8, color="black")+
  geom_scattermore(mapping=aes(UMAP_1,UMAP_2), pointsize = 7, color="white")+
  geom_scattermore(mapping=aes(UMAP_1,UMAP_2, color=annotation_level_4), pointsize=4)

p <- p+
    scale_colour_manual(values=col_cells)+
    ylab("UMAP2")+xlab("UMAP1")+
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 5),
          legend.title = element_text(colour="white", size=3),
          panel.background = element_rect(colour = "black", size=0.5),
          axis.text.x = element_text(colour="black"), 
          axis.text.y = element_text(colour="black"))+
    coord_fixed()+
  ggtitle("Reference dataset")+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  theme(legend.key.size = unit(0.5, "cm"))
p


## Get plot of BDNF:

# Plot the single cell data
f <- FeaturePlot(scRef, feature="SPP1", order=T)


col=colorRampPalette(c("#FFFFFF", RColorBrewer::brewer.pal(9, "Reds")))
library(scattermore)

p <- 
  ggplot(data=f$data %>% arrange((SPP1)))+
  geom_scattermore(mapping=aes(UMAP_1,UMAP_2), pointsize = 8, color="black")+
  geom_scattermore(mapping=aes(UMAP_1,UMAP_2), pointsize = 7, color="white")+
  geom_scattermore(mapping=aes(UMAP_1,UMAP_2, color=SPP1), pointsize=5, pixels = c(1000, 1000))+
  geom_point(mapping=aes(UMAP_1,UMAP_2, color=SPP1), size=0.5)

p <- p+
  scale_color_gradientn(colours = col(50), name="", 
                       limits=c(0,2), 
                       oob=scales::squish)+
    guides(color = guide_colourbar(barwidth = 0.3, barheight = 8, ticks =F, frame.colour="black"), label=F)+
    ylab("UMAP2")+xlab("UMAP1")+
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 5),
          legend.title = element_text(colour="white", size=3),
          panel.background = element_rect(colour = "black", size=0.5),
          axis.text.x = element_text(colour="black"), 
          axis.text.y = element_text(colour="black"))+
    coord_fixed()
p


## Create a module score:

# Create module Score

#scRef <- Seurat::AddModuleScore(scRef, features = list_genes, name=names(list_genes), assay = "SCT")
#names(list_genes)
library(dplyr)

modules <- seurat.obj@misc$vis$wgcna_modules$color %>% unique()
modules <- modules[2:length(modules)]
list_genes_imp <- map(.x = modules[1:length(modules)], .f = function(i) {
  importance <- 
    seurat.obj@misc$vis$wgcna_modules %>% 
    filter(module == i) %>% 
    dplyr::select(paste0("kME_", i))  # 显式调用 dplyr::select
  return(importance)
})


# 
# list_genes_imp <- map(.x=modules[1:length(modules)], .f=function(i){
#   
#   importance <- 
#     seurat.obj@misc$vis$wgcna_modules %>% 
#     filter(module==i) %>% 
#     select(paste0("kME_",i))
# 
#   return(importance)
#   
# })
names(list_genes_imp) <- modules

# get the amount of genes from each signature that is expressed

Module_signature <- map(list_genes_imp, .f=function(genes){
  
  # weighted mean expression
  mat <- Seurat::GetAssayData(scRef, assay = "RNA")
  mat <- mat[na.omit(match(rownames(genes), rownames(mat))),] %>% as.data.frame()
  
  mat <- mat*genes[rownames(mat), ]
  
  mean <- colMeans(mat)
  out <- apply(mat, 2, function(x){ length(which(x!=0))/length(x) })*100
  names(out) <- colnames(mat)
  mean=mean*out

  
  
  return(list(mean, out))
  
  
  
  
}, .progress = T)
names(Module_signature) <- paste0(names(list_genes_imp), 1:length(names(list_genes_imp)))

for(n in names(Module_signature)){
  scRef@meta.data[,n] <- Module_signature[[n]][[1]]
}

names(Module_signature)



## Plot modules 


module <- "turquoise2"

f <- FeaturePlot(scRef, features = module, order=T)

library(scattermore)
p <- ggplot(data=f$data)+
  geom_scattermore(mapping=aes(UMAP_1,UMAP_2), pointsize = 8, color="black")+
  geom_scattermore(mapping=aes(UMAP_1,UMAP_2), pointsize = 7, color="white")+
  geom_scattermore(mapping=aes(UMAP_1,UMAP_2, color=!!sym(module)), pointsize=2)

col=colorRampPalette(c("#FFFFFF", RColorBrewer::brewer.pal(9, "Greens")))
p <- p+
    scale_colour_gradientn(colours = col(50),
                           limit=c(0,5), 
                           oob = scales::squish, na.value = "white")+
    ylab("UMAP2")+xlab("UMAP1")+
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", size=0.5),
          axis.text.x = element_text(colour="black"), 
          axis.text.y = element_text(colour="black"))+
    coord_fixed()
p

##green3产生了NANANA




Module_score <- names(Module_signature)
Module_score <- Module_score[!Module_score %in% "green3"]

data_enrichment <- scRef@meta.data[,c("annotation_level_4", Module_score)]

enrichment <- 
  data_enrichment %>% 
  group_by(annotation_level_4) %>% 
  summarise_all(.funs = base::mean)

size <- map_dfc(1:length(Module_signature), ~ as.data.frame(Module_signature[[.x]][[2]]))
names(size) <- names(Module_signature)
size<-size[,-3]
size <- 
  size %>% 
  mutate(annotation_level_4=scRef$annotation_level_4) %>% 
  group_by(annotation_level_4) %>% 
  summarise_all(.funs = mean)


#tumors <- colors[colors$annotation_level_2 %in% c("Differentiated-like", "Stem-like"), ]$annotation_level_4
tumors <- colors$annotation_level_4
enrichment <- enrichment %>% filter(annotation_level_4 %in% tumors) %>% as.data.frame()
rownames(enrichment) <- enrichment$annotation_level_4
enrichment$annotation_level_4 <- NULL

size <- size %>% filter(annotation_level_4 %in% tumors) %>% as.data.frame()
rownames(size) <- size$annotation_level_4
size$annotation_level_4 <- NULL

enrichment_mat <- scales::rescale(as.matrix(enrichment), c(0,1)) %>% as.data.frame()
size <- scales::rescale(as.matrix(size), c(0,1)) %>% as.data.frame()
dim(enrichment)
dim(size)

enrichment <- 
  enrichment_mat %>%
  rownames_to_column(var = "gene") %>%
  gather(cell, count, -gene)

size <- 
  size %>%
  rownames_to_column(var = "gene") %>%
  gather(cell, count, -gene)

enrichment$size <- size$count

colors <- 
  colors %>% 
  filter(annotation_level_4!="Neuron")



col=colorRampPalette(c("#FFFFFF", RColorBrewer::brewer.pal(9, "Greens")))
#col=colorRampPalette(rev(RColorBrewer::brewer.pal(9, "PRGn")))
pheatmap::pheatmap(enrichment_mat %>% cor())


rank <- map_dbl(1:ncol(seurat.obj@misc$vis$MEs), ~cor(seurat.obj@misc$vis$MEs[.x], seurat.obj$HALLMARK.HYPOXIA))
names(rank) <- colnames(seurat.obj@misc$vis$MEs)
rank <- rank[modules]
names(rank) <- paste0(names(rank), 1:length(names(rank)))
rank <- rank[order(rank)]


library(ggcorrplot)
plot_1 <- 
  #ggcorrplot(enrichment, method = "circle", outline.color = "black")+
  ggplot(enrichment)+
  geom_point(aes(x = gene, y = cell, size = size, fill=count),colour="black",pch=21)+
  scale_y_discrete(limits = names(rank))+
  scale_x_discrete(limits = colors$annotation_level_4)+
  scale_fill_gradientn(colours = col(50), name="", limits=c(0,0.3), oob=scales::squish)+
  scale_size(range = c(0,3))+
  guides(fill = guide_colourbar(barwidth = 0.3, barheight = 8, ticks =F, frame.colour="black"), label=F)+
  ggtitle("Expression vs Metabolism")+
  xlab("")+ylab("")+
  theme_bw() +
    theme(panel.grid.major = e飞lement_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", size=0.5),
          axis.text.x = element_text(colour="black",
                                     angle = 90, vjust = 0.5, hjust=1), 
          axis.text.y = element_text(colour="black"))+
  coord_fixed()

plot_1
colnames(scRef@meta.data)
table(scRef$iCNV)
DimPlot(scRef, group.by  = "iCNV",raster=FALSE,cols =mycol)
DimPlot(scRef, group.by  = "annotation_level_3",raster=FALSE,cols =mycol)
table(scRef$annotation_level_3)
FeaturePlot(scRef, feature="SPP1", order=T)

# 方法2: 小提琴图（比较表达分布）
VlnPlot(scRef, features = "SPP1", group.by = "annotation_level_3", 
        cols = mycol, pt.size = 0.1, log = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 方法3: 点图（表达水平+表达比例）
DotPlot(scRef, features = "SPP1", group.by = "annotation_level_3", 
        cols = c("lightblue", "darkblue")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 方法4: 山脊图（密度分布）
RidgePlot(scRef, features = c("SPP1","CD44") ,group.by = "annotation_level_3", 
          cols = mycol, log = TRUE)
RidgePlot(scRef, features = "SPP1" ,group.by = "annotation_level_3", 
          cols = mycol, log = TRUE)
RidgePlot(scRef, features = "SPP1" ,group.by = "annotation_level_2", 
          cols = mycol, log = TRUE)
# 额外：获取具体的表达统计数据# 额外：获取具体的表达统计数据
# 计算每个细胞类型的平均表达量
avg_expression <- AverageExpression(scRef, features = "SPP1", 
                                    group.by = "annotation_level_3")
print(avg_expression)



