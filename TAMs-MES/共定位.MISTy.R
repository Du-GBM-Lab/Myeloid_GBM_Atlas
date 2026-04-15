
# MISTy
library(mistyR)
library(future)
#Seurat
library(Seurat)
library(SeuratObject)
# Data manipulation
library(tidyverse)
# Distances
library(distances)
library(qs)


gbm_UKF269<-qread('gbmUKF269.qs')


library(spacexr)
library(tidyverse)
options(future.globals.maxSize = 1000 * 1024^5)
options(stringsAsFactors = FALSE)
set.seed(123)
#读入reference并走标准流程 #细胞参考数据库构建 # 整体
scRNA <- readRDS('D:/ProgramData/Spatial Gene Expression/GBmap/data/azimuth_core_GBmap.rds')
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = VariableFeatures(object = .), verbose = FALSE) %>%
  RunPCA(features = VariableFeatures(object = .), verbose = FALSE)
sc_counts <- as.matrix(scRNA[['RNA']]@counts)
sc_nUMI = colSums(sc_counts)
table(scRNA$annotation_level_3)
scRNA$celltype<-scRNA$annotation_level_3
scRNA$celltype<-gsub('CD4/CD8','CD4-CD8',scRNA$celltype)
table(scRNA$celltype)
cellType=data.frame(barcode=colnames(scRNA), celltype=scRNA$celltype)
names(cellType) = c('barcode', 'cell_type')
cell_types = cellType$cell_type; names(cell_types) <- cellType$barcode 
cell_types =as.factor(cell_types)
reference = Reference(sc_counts, cell_types, sc_nUMI,min_UMI = 50)

spatial_count <- gbm_UKF269@assays$Spatial@counts
head(spatial_count)
#获取坐标
coords <- gbm_UKF269@images[["UKF269"]]@coordinates
head(coords)
#只保留row,col列
coords<-coords[,2:3]

#创建query 
query  <- SpatialRNA(coords, spatial_count)

#运行RCTD开始反卷积
RCTD <- create.RCTD(query, reference, max_cores = 1,CELL_MIN_INSTANCE = 20)

#10X Visium等低分辨率使用full模式(高分辨率可使用doublet)
RCTD <- run.RCTD(RCTD, doublet_mode = "full") 
qsave(RCTD,file = 'RCTD.qs')
RCTD<-qread('RCTD.qs')
#标准化细胞占比权重
weights <- RCTD@results$weights
norm_weights <- normalize_weights(weights)


# Add RCTD results to Seurat object
gbm_UKF269 <- AddMetaData(gbm_UKF269, metadata = RCTD@results$weights)

# Normalize weights
weights <- RCTD@results$weights
norm_weights <- normalize_weights(weights)

# Add RCTD results as a new assay
gbm_UKF269[["rctd_full"]] <- CreateAssayObject(data = t(as.matrix(norm_weights)))
if (length(gbm_UKF269@assays$rctd_full@key) == 0) {
  gbm_UKF269@assays$rctd_full@key <- "rctd_full_"
}

# Plotting
DefaultAssay(gbm_UKF269) <- "rctd_full"
colnames(gbm_UKF269@meta.data)
cell_types <- c('MES.like','AC.like','NPC.like','OPC.like','TAM.BDM','TAM.MG','Mono','DC')

SpatialFeaturePlot(gbm_UKF269, features = cell_types, pt.size.factor = 1.6, 
                   ncol = 4, crop = TRUE)

# for (img in image_names) {
#   plot <- SpatialFeaturePlot(gbm_UKF269, features = cell_types, pt.size.factor = 1.6, 
#                              ncol = 3, crop = TRUE, images = img)
#   ggsave(paste0("plots/", img, "_rctd_full.png"), plot, width = 12, height = 8, dpi = 300)
# }
library(pheatmap)
library(progeny)
# PROGENy analysis
gbm_UKF269_progeny <- progeny(gbm_UKF269, scale=FALSE, organism="Human", top=1000, 
                        perm=1, return_assay = TRUE, assay="Spatial")
gbm_UKF269_progeny <- ScaleData(gbm_UKF269_progeny, assay = "progeny")

progeny_scores_df <- as.data.frame(t(GetAssayData(gbm_UKF269_progeny, slot = "scale.data", 
                                                  assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell)

CellsClusters <- data.frame(Cell = names(Idents(gbm_UKF269)), 
                            CellType = as.character(Idents(gbm_UKF269)), 
                            stringsAsFactors = FALSE)

progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

summarized_progeny_scores <- progeny_scores_df %>% 
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))

summarized_progeny_scores_df <- summarized_progeny_scores %>%
  select(-std) %>%   
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 

# Heatmap of PROGENy scores
paletteLength <- 100
myColor <- colorRampPalette(c("#008080", "white","#FFA500"))(paletteLength)

progenyBreaks <- c(seq(min(summarized_progeny_scores_df), 0, 
                       length.out=ceiling(paletteLength/2) + 1),
                   seq(max(summarized_progeny_scores_df)/paletteLength, 
                       max(summarized_progeny_scores_df), 
                       length.out=floor(paletteLength/2)))

progeny_hmap <- pheatmap(t(summarized_progeny_scores_df[,-1]), fontsize=14, 
                         fontsize_row = 10, color=myColor, breaks = progenyBreaks, 
                         main = "PROGENy (1000)", angle_col = 45,
                         treeheight_col = 0,  border_color = NA)

# Save the heatmap
ggsave("progeny_heatmap.pdf", progeny_hmap, width = 8, height = 5, dpi = 700)






# 1. 提取RCTD反卷积结果
# -----------------------------------------------------------------
# 获取归一化的细胞类型比例矩阵 (spots × cell types)
celltype_composition <- as.data.frame(RCTD@results$weights)
colnames(celltype_composition) <- gsub(" ", "_", colnames(celltype_composition), fixed = TRUE)
colnames(celltype_composition) <- gsub("-", "_", colnames(celltype_composition), fixed = TRUE)
rownames(celltype_composition) <- colnames(gbm_UKF269)  # 确保spot名称一致

# 2. 准备空间坐标
# -----------------------------------------------------------------
coords <- GetTissueCoordinates(gbm_UKF269, 
                               cols = c("imagerow", "imagecol"), 
                               scale = NULL) %>%
  select(row = imagerow, col = imagecol)  # 重命名列以符合mistyR要求

# 3. 计算paraview最优半径
# -----------------------------------------------------------------
geom_dist <- as.matrix(distances(coords))          # 计算所有spot间欧氏距离
dist_nn <- apply(geom_dist, 1, function(x) sort(x)[2])  # 获取每个spot到最近邻的距离
paraview_radius <- ceiling(mean(dist_nn) + sd(dist_nn)) # 均值+标准差作为半径
print(paste("Paraview radius:", paraview_radius))  # 检查半径值（通常50-100μm）

# 4. 构建MISTy视图
# -----------------------------------------------------------------
misty_views <- create_initial_view(celltype_composition) %>%
  add_paraview(
    positions = coords,
    l = paraview_radius,        # 上一步计算的半径
    family = "gaussian")

# 5. 运行MISTy模型
# -----------------------------------------------------------------
run_misty(
  misty_views,
  results.folder = "./misty_results",  # 结果存储路径
  cached = FALSE                       # 强制重新计算
)




# 6. 收集与可视化结果
# -----------------------------------------------------------------
misty_results <- collect_results("./misty_results")

misty_results %>%
  plot_improvement_stats("multi.R2") %>% 
  plot_improvement_stats("gain.R2")

misty_results %>%plot_interaction_heatmap(view = "intra", clean = F)
misty_results %>%plot_interaction_heatmap(view = "intra", clean = F)


misty_results$importances.aggregated %>%
  filter(view == "intra", Predictor == "TAM_BDM") %>%
  arrange(-Importance)
head(misty_results)
misty_results$importances.aggregated %>%
  filter(view == "intra", Predictor == "TAM_MG") %>%
  arrange(-Importance)
head(misty_results)
misty_results$importances.aggregated %>%
  filter(view == "intra", Predictor == "Mono") %>%
  arrange(-Importance)
head(misty_results)
misty_results$importances.aggregated %>%
  filter(view == "intra", Predictor == "DC") %>%
  arrange(-Importance)
head(misty_results)
SpatialFeaturePlot(gbm_UKF269, keep.scale = NULL, features = c("MES.like" ,"TAM.BDM",'TAM.MG','OPC.like','Mono','Oligodendrocyte'), image.alpha = 0)

paraview_radius
misty_results %>% plot_interaction_heatmap(view = "para.70", clean = TRUE, 
                                           trim = 1.75, trim.measure = "gain.R2",
                                           cutoff = 0.1) 


SpatialFeaturePlot(gbm_UKF269, keep.scale = NULL, features = c("MES-like" ,"TAM-BDM"), image.alpha = 0)



misty_results %>% plot_view_contributions()
misty_results %>% plot_interaction_heatmap(view = "intra", cutoff = 0.8)
misty_results %>% plot_interaction_heatmap(view = "para.70", cutoff = 0.5)
misty_results %>% plot_contrast_heatmap("intra", "para.70", cutoff = 0.5)

misty_results %>% plot_interaction_communities("intra")
misty_results %>% plot_interaction_communities("para.70", cutoff = 0.5)



qsave(gbm_UKF269,file = 'gbm_UKF269_ISCHIA.qs')
