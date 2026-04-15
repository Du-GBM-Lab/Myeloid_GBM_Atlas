# Loading required packages
library(ISCHIA)
library(robustbase)
library(data.table)
library(ggplot2)
library(Seurat)
library(dplyr)
library(factoextra)
library(cluster)
library(showtext)
library(gridExtra)
library(pdftools)
library(qs)
library(dittoSeq)
library(spacexr)

# Set random seed for reproducibility
set.seed(123)

# Load data
gbm_UKF269<-qread('gbmUKF269.qs')
RCTD<-qread('RCTD.qs')
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
SpatialFeaturePlot(gbm_UKF269, features = cell_types, pt.size.factor = 1.6, 
                   ncol = 4, crop = TRUE,image.alpha = 0)

assay_matrix <- gbm_UKF269[["rctd_full"]]@data
norm_weights <- as.data.frame(t(assay_matrix))

# Elbow Method
k.values <- 1:20
wss_values <- sapply(k.values, function(k) kmeans(norm_weights, k, nstart = 10)$tot.withinss)

pdf("1_elbow_plot.pdf")
plot(k.values, wss_values, type = "b", pch = 19, frame = FALSE, 
     xlab = "Number of clusters K", ylab = "Total within-cluster sum of squares",
     main = "Elbow Method for Optimal K")
dev.off()

# Gap Statistic
gap_stat <- function(k) {
  km.res <- kmeans(norm_weights, k, nstart = 10)
  if (k == 1) return(NA)
  obs_disp <- sum(km.res$withinss)
  reference_disp <- mean(replicate(10, {
    km.null <- kmeans(matrix(rnorm(nrow(norm_weights) * ncol(norm_weights)), 
                             ncol = ncol(norm_weights)), k, nstart = 10)
    sum(km.null$withinss)
  }))
  log(reference_disp) - log(obs_disp)
}

gap_stat_values <- sapply(k.values, gap_stat)

pdf("2_gap_statistic_plot.pdf")
plot(k.values, gap_stat_values, type = "b", pch = 19, frame = FALSE, 
     xlab = "Number of Clusters (K)", ylab = "Gap Statistic",
     main = "Gap Statistic: Determining Optimal K")
dev.off()

# Calinski-Harabasz Index
calinski_harabasz_index <- function(data, labels) {
  num_clusters <- length(unique(labels))
  num_points <- nrow(data)
  
  # 使用aggregate计算每个簇的中心点，更稳健地处理维度
  # .groups参数确保输出是数据框
  centroids <- aggregate(data, by = list(cluster = labels), FUN = mean)
  # 移除分组列，只保留数值中心点
  centroid_matrix <- as.matrix(centroids[, -1, drop = FALSE])
  rownames(centroid_matrix) <- centroids$cluster # 可选，设置行名便于理解
  
  # 计算全局中心点
  global_center <- colMeans(data)
  
  # 计算簇间离散度 (Between-cluster dispersion)
  between_disp <- sum(sapply(1:num_clusters, function(i) {
    cluster_size <- sum(labels == i)
    # 确保使用矩阵索引
    cluster_center <- centroid_matrix[i, , drop = FALSE]
    cluster_size * sum((cluster_center - global_center)^2)
  }))
  
  # 计算簇内离散度 (Within-cluster dispersion)
  within_disp <- sum(sapply(1:num_clusters, function(i) {
    cluster_points <- data[labels == i, , drop = FALSE]
    # 确保使用矩阵索引
    cluster_center <- centroid_matrix[i, , drop = FALSE]
    # 计算簇内所有点到其簇中心的距离平方和
    sum(apply(cluster_points, 1, function(point) sum((point - cluster_center)^2)))
  }))
  
  # 计算Calinski-Harabasz指数
  if (num_clusters == 1) {
    return(NA) # 当只有一个簇时，指数未定义
  } else {
    ch_index <- (between_disp / (num_clusters - 1)) / (within_disp / (num_points - num_clusters))
    return(ch_index)
  }
}

ch_values <- sapply(k.values, function(k) {
  km.res <- kmeans(norm_weights, k, nstart = 10)
  calinski_harabasz_index(norm_weights, km.res$cluster)
})

pdf("3_calinski_harabasz_plot.pdf")
plot(k.values, ch_values, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of Clusters (K)", ylab = "Calinski-Harabasz Index",
     main = "Calinski-Harabasz Index: Determining Optimal K")
dev.off()

# ISCHIA Analysis
pdf("4_composition_cluster_k_plot.pdf")
Composition.cluster.k(norm_weights, 20)
dev.off()

gbm_UKF269 <- Composition.cluster(gbm_UKF269, norm_weights, 12)
gbm_UKF269$cc_12 <- gbm_UKF269$CompositionCluster_CC

# Spatial Dimension Plot
image_names <- c("IU_PDA_T1", "IU_PDA_T2", "IU_PDA_HM2", "IU_PDA_HM2_2", "IU_PDA_NP2", 
                 "IU_PDA_T3", "IU_PDA_HM3", "IU_PDA_T4", "IU_PDA_HM4", "IU_PDA_HM5", 
                 "IU_PDA_T6", "IU_PDA_HM6", "IU_PDA_LNM6", "IU_PDA_LNM7", "IU_PDA_T8", 
                 "IU_PDA_HM8", "IU_PDA_LNM8", "IU_PDA_T9", "IU_PDA_HM9", "IU_PDA_T10", 
                 "IU_PDA_HM10", "IU_PDA_LNM10", "IU_PDA_NP10", "IU_PDA_T11", "IU_PDA_HM11", 
                 "IU_PDA_NP11", "IU_PDA_T12", "IU_PDA_HM12", "IU_PDA_LNM12", "IU_PDA_HM13")

paletteMartin <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', 
                   '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff')

all_ccs <- unique(gbm_UKF269$CompositionCluster_CC)
color_mapping <- setNames(paletteMartin[1:length(all_ccs)], all_ccs)

pdf("5_spatial_plots_K12.pdf", width = 10, height = 7)

plot <- SpatialDimPlot(gbm_UKF269, group.by = "CompositionCluster_CC") +
    scale_fill_manual(values = color_mapping) +
    theme_minimal()

print(plot)

dev.off()

# Enriched Cell Types
save_cc_plot <- function(cc) {
  plot <- Composition_cluster_enrichedCelltypes(gbm_UKF269, cc, as.matrix(norm_weights))
  pdf_name <- paste0(cc, ".pdf")
  pdf(file = pdf_name)
  print(plot)
  dev.off()
}

ccs <- paste0("CC", 1:12)
for (cc in ccs) {
  save_cc_plot(cc)
}

pdf_files <- paste0("CC", 1:12, ".pdf")
pdf_combine(pdf_files, output = "6_enrichedCelltypes_CC_12.pdf")

# UMAP
gbm_UKF269.umap <- Composition_cluster_umap(gbm_UKF269, norm_weights)
pdf("7_umap_pie_chart.pdf")
print(gbm_UKF269.umap$umap.deconv.gg)
dev.off()

# Add UMAP to Seurat object
emb.umap <- gbm_UKF269.umap$umap.table
emb.umap$CompositionCluster_CC <- NULL
emb.umap$Slide <- NULL
emb.umap <- as.matrix(emb.umap)
colnames(emb.umap) <- c("UMAP1", "UMAP2")

gbm_UKF269[['umap.ischia12']] <- CreateDimReducObject(embeddings = emb.umap, key = 'umap.ischia12_', assay = 'rctd_full')

pdf("8_seurat_ischia_umap_12.pdf")
DimPlot(gbm_UKF269, reduction = "umap.ischia12", label = FALSE, group.by="cc_12")
dev.off()
colnames(gbm_UKF269@meta.data)
# Bar plots
pdf("9_barplot_SampVsorig_12.pdf", height=12, width=20)
dittoBarPlot(gbm_UKF269, "Metaprograms", group.by = "cc_12")
dev.off()

pdf("10_barplot_origVsSamp_12.pdf", height=10, width=10)
dittoBarPlot(gbm_UKF269, "cc_12", group.by = "Metaprograms")
dev.off()

# Cell type co-occurrence
CC4.celltype.cooccur <- spatial.celltype.cooccurence(spatial.object=gbm_UKF269, deconv.prob.mat=norm_weights, 
                                                     COI="CC4", prob.th= 0.05, 
                                                     Condition=unique(gbm_UKF269$orig.ident))
pdf("11_celltype_cooccurrence_CC4.pdf")
plot.celltype.cooccurence(CC4.celltype.cooccur)
dev.off()
CC10.celltype.cooccur <- spatial.celltype.cooccurence(spatial.object=gbm_UKF269, deconv.prob.mat=norm_weights, 
                                                     COI="CC10", prob.th= 0.05, 
                                                     Condition=unique(gbm_UKF269$orig.ident))
pdf("11_celltype_cooccurrence_CC4.pdf")
plot.celltype.cooccurence(CC10.celltype.cooccur)
dev.off()






