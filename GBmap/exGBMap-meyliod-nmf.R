.libPaths(c('~/R/x86_64-pc-linux-gnu-library/4.4',
            '/refdir/Rlib',
            '/usr/local/lib/R/library'))
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.4/")
getwd()
.libPaths()
library(Seurat)
library(patchwork)
library(ggplot2)
library(bigmemory)
library(patchwork)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(ggsci)
library(ggrastr)
library(clusterProfiler)
library(GeneNMF)
library(BiocParallel)
library(viridis)
library(qs)
custom_magma <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))

bpparam()  # 检查当前的并行参数
plan("multisession", workers = 8)
options(future.globals.maxSize = 100000000 * 1024^5)
options(stringsAsFactors = FALSE)
set.seed(123)
mycol <- c(pal_d3()(7),pal_aaas()(7),pal_uchicago()(7),pal_jama()(7))
meyloid_altlas <- qread("meyloid_stand.qs") 

colnames(meyloid_altlas@meta.data)
table(meyloid_altlas$annotation_level_1)
table(meyloid_altlas$annotation_level_2)
table(meyloid_altlas$annotation_level_3)
table(meyloid_altlas$annotation_level_4)


DimPlot(meyloid_altlas,group.by = 'annotation_level_2')

# table(meyloid_altlas$patient)
# DimPlot(meyloid_altlas,group.by = 'patient')
# select_cell <- rownames(meyloid_altlas@meta.data)[as.vector(meyloid_altlas@meta.data[,"annotation_level_1"]) %in% c("Neoplastic")]
# seu <- subset(meyloid_altlas, cells = select_cell)
seu<-meyloid_altlas
rm(meyloid_altlas)
colnames(seu@meta.data)

DimPlot(seu, group.by="annotation_level_3")
DimPlot(seu, group.by="annotation_level_4")
DimPlot(seu,reduction = 'tsne')
DefaultAssay(seu) <- "RNA"
colnames(seu@meta.data)
seu.list <- SplitObject(meyloid_altlas, split.by = "donor_id")

geneNMF.programs <- multiNMF(seu.list,k=4:9)
geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs,
                                        metric = "cosine",
                                        weight.explained = 0.5,
                                        nMP=15) #5cell 5common 4mp
saveRDS(geneNMF.programs,file = 'geneNMF.programs.rds')
saveRDS(geneNMF.metaprograms,file = "geneNMF.metaprograms.rds")
#geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs, nMP=4)
ph <- plotMetaPrograms(geneNMF.metaprograms,palette = custom_magma,
                       similarity.cutoff = c(0.1,1))

ph
##plotMetaPrograms-meyloid_exmeyloid_altlas-stand

plotMetaPrograms(geneNMF.metaprograms,similarity.cutoff = c(0.1,1),
                 palette = c('white','#bfd3e6','#9ebcda','#8c96c6',
                             '#8c6bb1','#88419d','#810f7c','#4d004b'))

geneNMF.metaprograms_jaccard <- getMetaPrograms(geneNMF.programs,
                                                metric = "jaccard",
                                                weight.explained = 0.5,nMP=15) #14出13个
ph_j <- plotMetaPrograms(geneNMF.metaprograms_jaccard,palette = custom_magma,similarity.cutoff = c(0, 0.3))
ph_j
dev.off()
plotMetaPrograms(geneNMF.metaprograms_jaccard,palette = custom_magma,similarity.cutoff = c(0, 0.1))
plotMetaPrograms(geneNMF.metaprograms_jaccard,palette = custom_magma,similarity.cutoff = c(0, 0.00001))

geneNMF.metaprograms_jaccard <- getMetaPrograms(geneNMF.programs,
                                                metric = "jaccard",
                                                weight.explained = 0.85,nMP=15) #14出13个
pdf("custom_plot.pdf", family = "Times", bg = "lightgray",width =25, height=16)  # 设置字体和背景色
plotMetaPrograms(geneNMF.metaprograms_jaccard,palette = custom_magma,similarity.cutoff = c(0, 0.1))
dev.off()  # 

geneNMF.metaprograms$metaprograms.metrics
# geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs,
#                                         metric = "cosine",
#                                         weight.explained = 0.5,
#                                         nMP=5,
#                                         min.confidence = 0.8)
# geneNMF.metaprograms$metaprograms.metrics
#   
# ph <- plotMetaPrograms(geneNMF.metaprograms,
#                        similarity.cutoff = c(0.1,1),
#                        scale = "none",
#                        palette = viridis(100, option = "A", direction = -1),
#                        annotation_colors = NULL,
#                        main = "Clustered Heatmap",
#                        show_rownames = FALSE,
#                        show_colnames = FALSE)
# ph
#rm(seu,seu.list)
geneNMF.metaprograms$metaprograms.metrics
lapply(geneNMF.metaprograms$metaprograms.genes, head)
geneNMF.metaprograms$metaprograms.genes.weights$MP1





library(msigdbr)
library(fgsea)
top_p <- lapply(geneNMF.metaprograms$metaprograms.genes, function(program) {
  runGSEA(program, universe=rownames(seu), category = "C5", subcategory = "GO:BP")
})
head(top_p$MP1)
head(top_p$MP2)
head(top_p$MP3)
head(top_p$MP4)

library(UCell)
mp.genes <- geneNMF.metaprograms$metaprograms.genes
seu <- AddModuleScore_UCell(seu, features = mp.genes, slot = "data",name = "", ncores=20)
saveRDS(seu,file = "AddModuleScore_UCell-seu.rds")
VlnPlot(seu, features=names(mp.genes), group.by = "annotation_level_2",
        pt.size = 0, ncol=3,raster=FALSE)
matrix <- seu@meta.data[,names(mp.genes)] 

#dimred <- scale(matrix)
dimred <- as.matrix(matrix)

colnames(dimred) <- paste0("MP_",seq(1, ncol(dimred)))
#New dim reduction
seu@reductions[["MPsignatures"]] <- new("DimReduc",
                                        cell.embeddings = dimred,
                                        assay.used = "RNA",
                                        key = "MP_",
                                        global = FALSE)
set.seed(123)
seu <- RunUMAP(seu, reduction="MPsignatures", dims=1:length(seu@reductions[["MPsignatures"]]),
               metric = "euclidean", reduction.name = "umap_MP")
colnames(seu@meta.data)
DimPlot(seu, reduction = "umap_MP", group.by = "patient") + theme(aspect.ratio = 1)
DimPlot(seu, reduction = "umap_MP", group.by = "annotation_level_2") + theme(aspect.ratio = 1)
DimPlot(seu, reduction = "umap_MP", group.by = "annotation_level_3") + theme(aspect.ratio = 1)
DimPlot(seu, reduction = "umap_MP", group.by = "celltype_original") + theme(aspect.ratio = 1)
DimPlot(seu, reduction = "umap_MP", group.by = "annotation_level_1") + theme(aspect.ratio = 1)

FeaturePlot(seu,features= "Systemic1")
FeaturePlot(seu,features= "MG_Inflamm2")
FeaturePlot(seu,features= "Complement3")
FeaturePlot(seu,features= "Scavenger4")
FeaturePlot(meyloid_altlas,features= "Systemic1")
FeaturePlot(meyloid_altlas,features= "MG_Inflamm2")
FeaturePlot(meyloid_altlas,features= "Complement3")
FeaturePlot(meyloid_altlas,features= "Scavenger4")

library(viridis)
FeaturePlot(seu, features = names(mp.genes), reduction = "umap_MP", ncol=3) &
  scale_color_viridis(option="A") &
  theme(aspect.ratio = 1, axis.text=element_blank(), axis.ticks=element_blank())




