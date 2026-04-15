library(Seurat)
library(SeuratData)
library(SeuratObject)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(forcats)
library(irGSEA)
library(msigdbr)
library(SeuratData)
library(RcppML)
library(tidyverse) 
library(infercnv) 
library(harmony) 
library(export) 
library(ggpubr)
plan("multisession", workers = 20)
options(future.globals.maxSize = 1000 * 1024^6)
options(stringsAsFactors = FALSE)

meyloid_SeuratObj<-readRDS('meyloid_SeuratObj.RDS')

table(meyloid_SeuratObj$stage)
seurat_obj = meyloid_SeuratObj[,meyloid_SeuratObj$stage %in% c("Normal","Early")]
meta <- seurat_obj@meta.data


###CNV_Early
infer_CNV_obj<-readRDS('./CNV_Early/run.final.infercnv_obj')
expr<-infer_CNV_obj@expr.data
expr[1:4,1:4]
data_cnv<-as.data.frame(expr)
dim(expr)
colnames(data_cnv)
rownames(data_cnv)


if(T){
  tmp1 = expr[,infer_CNV_obj@reference_grouped_cell_indices$`ref-1`]
  tmp2 = expr[,infer_CNV_obj@reference_grouped_cell_indices$`ref-2`]
  tmp= cbind(tmp1,tmp2)
  down=mean(rowMeans(tmp)) - 2 * mean( apply(tmp, 1, sd))
  up=mean(rowMeans(tmp)) + 2 * mean( apply(tmp, 1, sd))
  oneCopy=up-down
  oneCopy
  a1= down- 2*oneCopy
  a2= down- 1*oneCopy
  down;up
  a3= up +  1*oneCopy
  a4= up + 2*oneCopy 
  
  cnv_score_table<-infer_CNV_obj@expr.data
  cnv_score_table[1:4,1:4]
  cnv_score_mat <- as.matrix(cnv_score_table)
  
  # Scoring
  cnv_score_table[cnv_score_mat > 0 & cnv_score_mat < a2] <- "A" #complete loss. 2pts
  cnv_score_table[cnv_score_mat >= a2 & cnv_score_mat < down] <- "B" #loss of one copy. 1pts
  cnv_score_table[cnv_score_mat >= down & cnv_score_mat <  up ] <- "C" #Neutral. 0pts
  cnv_score_table[cnv_score_mat >= up  & cnv_score_mat <= a3] <- "D" #addition of one copy. 1pts
  cnv_score_table[cnv_score_mat > a3  & cnv_score_mat <= a4 ] <- "E" #addition of two copies. 2pts
  cnv_score_table[cnv_score_mat > a4] <- "F" #addition of more than two copies. 2pts
  
  # Check
  table(cnv_score_table[,1])
  # Replace with score 
  cnv_score_table_pts <- cnv_score_mat
  rm(cnv_score_mat)
  # 
  cnv_score_table_pts[cnv_score_table == "A"] <- 2
  cnv_score_table_pts[cnv_score_table == "B"] <- 1
  cnv_score_table_pts[cnv_score_table == "C"] <- 0
  cnv_score_table_pts[cnv_score_table == "D"] <- 1
  cnv_score_table_pts[cnv_score_table == "E"] <- 2
  cnv_score_table_pts[cnv_score_table == "F"] <- 2
  
  cnv_score_table_pts[1:4,1:4]
  str(  as.data.frame(cnv_score_table_pts[1:4,1:4])) 
  cell_scores_CNV <- as.data.frame(colSums(cnv_score_table_pts))
  
  colnames(cell_scores_CNV) <- "cnv_score" 
}

head(cell_scores_CNV) 
score=cell_scores_CNV
head(score)
meta$totalCNV = score[match(colnames(seurat_obj),
                            rownames(score)),1] 

ggplot(meta, aes(x=cell_type  , y=log10(log10(totalCNV)), fill=cell_type  )) +
  geom_boxplot()   
ggplot(meta, aes(x = cell_type, y = totalCNV, fill = cell_type)) +
  geom_boxplot(
    width = 0.6,               # 缩窄箱体宽度
    outlier.shape = NA,        # 隐藏异常值
    alpha = 0.8                # 增加透明度
  ) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.3) +  # 添加抖动点
  scale_y_continuous(limits = c(3600, 3640), breaks = seq(3600, 3640, by = 5)) +  # 精细刻度
  theme_minimal()



table(meyloid_SeuratObj$stage)
seurat_obj = meyloid_SeuratObj[,meyloid_SeuratObj$stage %in% c("Normal","Medium")]
meta <- seurat_obj@meta.data


###CNV_Medium
infer_CNV_obj<-readRDS('./CNV_Medium//run.final.infercnv_obj')
expr<-infer_CNV_obj@expr.data
expr[1:4,1:4]
data_cnv<-as.data.frame(expr)
dim(expr)
colnames(data_cnv)
rownames(data_cnv)


if(T){
  tmp1 = expr[,infer_CNV_obj@reference_grouped_cell_indices$`ref-1`]
  tmp2 = expr[,infer_CNV_obj@reference_grouped_cell_indices$`ref-2`]
  tmp= cbind(tmp1,tmp2)
  down=mean(rowMeans(tmp)) - 2 * mean( apply(tmp, 1, sd))
  up=mean(rowMeans(tmp)) + 2 * mean( apply(tmp, 1, sd))
  oneCopy=up-down
  oneCopy
  a1= down- 2*oneCopy
  a2= down- 1*oneCopy
  down;up
  a3= up +  1*oneCopy
  a4= up + 2*oneCopy 
  
  cnv_score_table<-infer_CNV_obj@expr.data
  cnv_score_table[1:4,1:4]
  cnv_score_mat <- as.matrix(cnv_score_table)
  
  # Scoring
  cnv_score_table[cnv_score_mat > 0 & cnv_score_mat < a2] <- "A" #complete loss. 2pts
  cnv_score_table[cnv_score_mat >= a2 & cnv_score_mat < down] <- "B" #loss of one copy. 1pts
  cnv_score_table[cnv_score_mat >= down & cnv_score_mat <  up ] <- "C" #Neutral. 0pts
  cnv_score_table[cnv_score_mat >= up  & cnv_score_mat <= a3] <- "D" #addition of one copy. 1pts
  cnv_score_table[cnv_score_mat > a3  & cnv_score_mat <= a4 ] <- "E" #addition of two copies. 2pts
  cnv_score_table[cnv_score_mat > a4] <- "F" #addition of more than two copies. 2pts
  
  # Check
  table(cnv_score_table[,1])
  # Replace with score 
  cnv_score_table_pts <- cnv_score_mat
  rm(cnv_score_mat)
  # 
  cnv_score_table_pts[cnv_score_table == "A"] <- 2
  cnv_score_table_pts[cnv_score_table == "B"] <- 1
  cnv_score_table_pts[cnv_score_table == "C"] <- 0
  cnv_score_table_pts[cnv_score_table == "D"] <- 1
  cnv_score_table_pts[cnv_score_table == "E"] <- 2
  cnv_score_table_pts[cnv_score_table == "F"] <- 2
  
  cnv_score_table_pts[1:4,1:4]
  str(  as.data.frame(cnv_score_table_pts[1:4,1:4])) 
  cell_scores_CNV <- as.data.frame(colSums(cnv_score_table_pts))
  
  colnames(cell_scores_CNV) <- "cnv_score" 
}

head(cell_scores_CNV) 
score=cell_scores_CNV
head(score)
meta$totalCNV = score[match(colnames(seurat_obj),
                            rownames(score)),1] 

ggplot(meta, aes(x=cell_type  , y=log10(log10(totalCNV)), fill=cell_type  )) +
  geom_boxplot()   
ggplot(meta, aes(x = cell_type, y = totalCNV, fill = cell_type)) +
  geom_boxplot(
    width = 0.6,               # 缩窄箱体宽度
    outlier.shape = NA,        # 隐藏异常值
    alpha = 0.8                # 增加透明度
  ) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.3) +  # 添加抖动点
  scale_y_continuous(limits = c(3600, 3640), breaks = seq(3600, 3640, by = 5)) +  # 精细刻度
  theme_minimal()






table(meyloid_SeuratObj$stage)
seurat_obj = meyloid_SeuratObj[,meyloid_SeuratObj$stage %in% c("Normal","Late")]
meta <- seurat_obj@meta.data


###CNV_Medium
infer_CNV_obj<-readRDS('./CNV_Late//run.final.infercnv_obj')
expr<-infer_CNV_obj@expr.data
expr[1:4,1:4]
data_cnv<-as.data.frame(expr)
dim(expr)
colnames(data_cnv)
rownames(data_cnv)


if(T){
  tmp1 = expr[,infer_CNV_obj@reference_grouped_cell_indices$`ref-1`]
  tmp2 = expr[,infer_CNV_obj@reference_grouped_cell_indices$`ref-2`]
  tmp= cbind(tmp1,tmp2)
  down=mean(rowMeans(tmp)) - 2 * mean( apply(tmp, 1, sd))
  up=mean(rowMeans(tmp)) + 2 * mean( apply(tmp, 1, sd))
  oneCopy=up-down
  oneCopy
  a1= down- 2*oneCopy
  a2= down- 1*oneCopy
  down;up
  a3= up +  1*oneCopy
  a4= up + 2*oneCopy 
  
  cnv_score_table<-infer_CNV_obj@expr.data
  cnv_score_table[1:4,1:4]
  cnv_score_mat <- as.matrix(cnv_score_table)
  
  # Scoring
  cnv_score_table[cnv_score_mat > 0 & cnv_score_mat < a2] <- "A" #complete loss. 2pts
  cnv_score_table[cnv_score_mat >= a2 & cnv_score_mat < down] <- "B" #loss of one copy. 1pts
  cnv_score_table[cnv_score_mat >= down & cnv_score_mat <  up ] <- "C" #Neutral. 0pts
  cnv_score_table[cnv_score_mat >= up  & cnv_score_mat <= a3] <- "D" #addition of one copy. 1pts
  cnv_score_table[cnv_score_mat > a3  & cnv_score_mat <= a4 ] <- "E" #addition of two copies. 2pts
  cnv_score_table[cnv_score_mat > a4] <- "F" #addition of more than two copies. 2pts
  
  # Check
  table(cnv_score_table[,1])
  # Replace with score 
  cnv_score_table_pts <- cnv_score_mat
  rm(cnv_score_mat)
  # 
  cnv_score_table_pts[cnv_score_table == "A"] <- 2
  cnv_score_table_pts[cnv_score_table == "B"] <- 1
  cnv_score_table_pts[cnv_score_table == "C"] <- 0
  cnv_score_table_pts[cnv_score_table == "D"] <- 1
  cnv_score_table_pts[cnv_score_table == "E"] <- 2
  cnv_score_table_pts[cnv_score_table == "F"] <- 2
  
  cnv_score_table_pts[1:4,1:4]
  str(  as.data.frame(cnv_score_table_pts[1:4,1:4])) 
  cell_scores_CNV <- as.data.frame(colSums(cnv_score_table_pts))
  
  colnames(cell_scores_CNV) <- "cnv_score" 
}

head(cell_scores_CNV) 
score=cell_scores_CNV
head(score)
meta$totalCNV = score[match(colnames(seurat_obj),
                            rownames(score)),1] 

ggplot(meta, aes(x=cell_type  , y=log10(log10(totalCNV)), fill=cell_type  )) +
  geom_boxplot()   
ggplot(meta, aes(x = cell_type, y = totalCNV, fill = cell_type)) +
  geom_boxplot(
    width = 0.6,               # 缩窄箱体宽度
    outlier.shape = NA,        # 隐藏异常值
    alpha = 0.8                # 增加透明度
  ) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.3) +  # 添加抖动点
  scale_y_continuous(limits = c(3600, 3640), breaks = seq(3600, 3640, by = 5)) +  # 精细刻度
  theme_minimal()


