
rm(list = ls())
library(tidyverse)
library(CellChat)
library(Seurat)
library(qs)
# Read data
visium.brain <- qread('D:/ProgramData/Glioma-immunelncRNA/TAMs-MES/gbmUKF269.qs')
colnames(visium.brain@meta.data)
Idents(visium.brain) <- "L3_second_type"
color.use <- scPalette(nlevels(visium.brain)); names(color.use) <- levels(visium.brain)
# Prepare input data for CelChat analysis
data.input = Seurat::GetAssayData(visium.brain, slot = "data", assay = "SCT") # normalized data matrix


meta = data.frame(labels = Seurat::Idents(visium.brain), samples = "sample1", row.names = names(Seurat::Idents(visium.brain))) # manually create a dataframe consisting of the cell labels
meta$samples <- factor(meta$samples)
unique(meta$labels) # check the cell labels
spatial.locs = Seurat::GetTissueCoordinates(visium.brain, scale = NULL, cols = c("imagerow", "imagecol"))

scalefactors <- jsonlite::fromJSON(
  txt = file.path(
    "D:/ProgramData/Glioma-immunelncRNA/TAMs-MES/UKF269_T_ST/outs/spatial",
    "scalefactors_json.json"
  )
)
spot.size = 65 # the theoretical spot size (um) in 10X Visium
conversion.factor = spot.size/scalefactors$spot_diameter_fullres
spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)

d.spatial <- computeCellDistance(coordinates = spatial.locs, ratio = spatial.factors$ratio, tol = spatial.factors$tol)
min(d.spatial[d.spatial!=0]) # this value should approximately equal 100um for 10X Visium data

#创建 CellChat 对象
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                           datatype = "spatial", coordinates = spatial.locs, spatial.factors = spatial.factors)
cellchat

#设置配体-受体相互作用数据库
CellChatDB <- CellChatDB.human # use CellChatDB.human if running on human data

showDatabaseCategory(CellChatDB)

dplyr::glimpse(CellChatDB$interaction)

# 使用CellChatDB的一个子集进行细胞-细胞通讯分析
CellChatDB.use <- subsetDB(CellChatDB,
                           search = "Secreted Signaling",
                           key = "annotation") # use Secreted Signaling

# 使用CellChatDB的所有子集进行分析
CellChatDB.use <- CellChatDB 

# 设定数据库所用数据集
cellchat@DB <- CellChatDB.use

#对子信号基因的表达数据进行子集提取，以节省计算成本
cellchat <- subsetData(cellchat) # 即便采用完整数据库，这一步骤仍然必不可少
future::plan("multisession", workers = 4) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat, variable.both = F)
#> The number of highly variable ligand-receptor pairs used for signaling inference is 1469 

#将基因表达数据投影到蛋白质–蛋白质相互作用网络（可选：在运行时，用户应在 computeCommunProb() 函数中将 raw.use = FALSE，以使用投影后的数据）
cellchat <- smoothData(cellchat, adj=PPI.human) # PPI.mouse

# 这一步骤非常缓慢
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1,
                              distance.use = TRUE, 
                              interaction.range = 250, 
                              scale.distance = 0.01,
                              contact.dependent = TRUE, 
                              contact.range = 100)

# 如果某些细胞组中的细胞数量较少，用户可以选择过滤掉细胞间的通信。默认情况下，每个细胞组中需要至少有 10 个细胞才能进行细胞间通信。
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
#可视化聚合的细胞间通讯网
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Blues")
#> Do heatmap based on a single object

netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Blues")
pathways.show <- c("SPP1") 
# Circle plot
par(mfrow=c(1,1), xpd = TRUE) # `xpd = TRUE` should be added to show the title
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial", edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 3.5)

#计算并可视化网络中心性得分
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # 槽 ‘netP’ 表示推断出的信号通路的细胞间通讯网络
# 使用热图可视化计算出的中心性分数，便于快速识别各细胞群在主要信号传导中的作用
par(mfrow=c(1,1))
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

# 用户在可视化信号网络时，可以在空间转录组图上展示这些信息，例如，圆圈越大表示接收到的信号越强
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial", edge.width.max = 2, alpha.image = 0.2, vertex.weight = "incoming", vertex.size.max = 4, vertex.label.cex = 3.5)

# 在组织中可视化基因表达分布
spatialFeaturePlot(cellchat, features = c("SPP1","CD44"), point.size = 1.2, color.heatmap = "Reds", direction = 1)

# 输入一个配体–受体对
spatialFeaturePlot(cellchat, pairLR.use = "IGF1_IGF1R", point.size = 0.5, do.binary = FALSE, cutoff = 0.05, enriched.only = F, color.heatmap = "Reds", direction = 1)

# 输入一个配体–受体对，并以二值形式显示其表达情况
spatialFeaturePlot(cellchat, pairLR.use = "IGF1_IGF1R", point.size = 1, do.binary = TRUE, cutoff = 0.05, enriched.only = F, color.heatmap = "Reds", direction = 1)
qsave(cellchat,file = 'gbm_UKF269_cellchat-L3_second_type.qs')
cellchat <- qread('D:/ProgramData/Glioma-immunelncRNA/TAMs-MES/gbm_UKF269_cellchat-L3_second_type.qs')



