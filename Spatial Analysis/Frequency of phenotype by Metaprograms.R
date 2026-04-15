# 加载必要包
library(Seurat)
library(pheatmap)
library(dplyr)
library(qs)

# select_cell <- rownames(stRNA234@meta.data)[as.vector(stRNA234@meta.data[,"L2_second_type"]) %in% c("Myeloid")]
# seu<- subset(stRNA234, cells = select_cell)
library(pheatmap)
library(tidyr)
library(tibble)
library(ggplot2)
stRNA234<-qread('stRNA-anno-level234_MPs.qs')
colnames(stRNA234@meta.data)
meta_df <- stRNA234@meta.data
# 计算Metaprograms
freq_table_4 <- as.data.frame(table(meta_df$layer, meta_df$Metaprograms))
colnames(freq_table_4) <- c("layer", "Metaprograms", "Count")

# 绘制堆叠条形图
ggplot(freq_table_4, aes(x =layer , y = Count, fill = Metaprograms)) +
  geom_bar(stat = "identity", position = "fill") +  # position="fill"转换为比例
  labs(y = "Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right") +
  scale_fill_viridis_d(option = "D")  # 使用颜色渐变
ggplot(freq_table_4, aes(x =Metaprograms , y = Count, fill = layer)) +
  geom_bar(stat = "identity", position = "fill") +  # position="fill"转换为比例
  labs(y = "Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right") +
  scale_fill_viridis_d(option = "D")  # 使用颜色渐变
 
# 标准化（按行或列） heatmap of layer by Metaprograms
heatmap_data_norm_4 <- prop.table(table(meta_df$layer, meta_df$Metaprograms), margin = 2)  # 按列标准化

# 绘制热图
p<-pheatmap(heatmap_data_norm_4,
            color = colorRampPalette(c("white", "red"))(100),
            cluster_rows = TRUE,
            cluster_cols = TRUE,
            #main = "Frequency of L4_second_type by mp",
            angle_col = '45',
            fontsize_row = 15,
            fontsize_col = 15)
ggsave('Frequency of layer by Metaprograms.pdf',plot = p,width = 6, height = 12, units = "in")

# 计算交叉频数表L4
freq_table_4 <- as.data.frame(table(meta_df$L4_second_type, meta_df$ivygap))
colnames(freq_table_4) <- c("L4_type", "ivygap", "Count")

# 绘制堆叠条形图
ggplot(freq_table_4, aes(x =L4_type , y = Count, fill = ivygap)) +
  geom_bar(stat = "identity", position = "fill") +  # position="fill"转换为比例
  labs(y = "Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right") +
  scale_fill_viridis_d(option = "D")  # 使用颜色渐变

# 标准化（按行或列）
heatmap_data_norm_4 <- prop.table(table(meta_df$L4_second_type, meta_df$ivygap), margin = 2)  # 按列标准化

# 绘制热图
p<-pheatmap(heatmap_data_norm_4,
         color = colorRampPalette(c("white", "red"))(100),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         #main = "Frequency of L4_second_type by mp",
         angle_col = '45',
         fontsize_row = 15,
         fontsize_col = 15)
ggsave('Frequency of L4_second_type by ivygap.pdf',plot = p,width = 6, height = 12, units = "in")
# 计算交叉频数表L3
freq_table_3 <- as.data.frame(table(meta_df$L3_second_type, meta_df$ivygap))
colnames(freq_table_3) <- c("L3_type", "ivygap", "Count")

# 绘制堆叠条形图
ggplot(freq_table_3, aes(x =L3_type , y = Count, fill = ivygap)) +
  geom_bar(stat = "identity", position = "fill") +  # position="fill"转换为比例
  labs(y = "Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right") +
  scale_fill_viridis_d(option = "D")  # 使用颜色渐变

# 标准化（按行或列）
heatmap_data_norm_3 <- prop.table(table(meta_df$L3_second_type, meta_df$ivygap), margin = 2)  # 按列标准化

# 绘制热图
pheatmap(heatmap_data_norm_3,
         color = colorRampPalette(c("white", "red"))(100),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         #main = "Frequency of L4_second_type by mp",
         angle_col = '45',
         fontsize_row = 15,
         fontsize_col = 15)




# 计算交叉频数表L2
freq_table_2 <- as.data.frame(table(meta_df$L2_second_type, meta_df$ivygap))
colnames(freq_table_2) <- c("L2_type", "ivygap", "Count")

# 绘制堆叠条形图
ggplot(freq_table_2, aes(x =L2_type , y = Count, fill = ivygap)) +
  geom_bar(stat = "identity", position = "fill") +  # position="fill"转换为比例
  labs(y = "Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right") +
  scale_fill_viridis_d(option = "D")  # 使用颜色渐变

# 标准化（按行或列）
heatmap_data_norm_2 <- prop.table(table(meta_df$L2_second_type, meta_df$ivygap), margin = 2)  # 按列标准化

# 绘制热图
pheatmap(heatmap_data_norm_2,
         color = colorRampPalette(c("white", "red"))(100),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         #main = "Frequency of L4_second_type by mp",
         angle_col = '45',
         fontsize_row = 15,
         fontsize_col = 15)

###########################
# 计算交叉频数表L4
freq_table_4 <- as.data.frame(table(meta_df$L4_second_type, meta_df$layer))
colnames(freq_table_4) <- c("L4_type", "layer", "Count")

# 绘制堆叠条形图
ggplot(freq_table_4, aes(x =L4_type , y = Count, fill = layer)) +
  geom_bar(stat = "identity", position = "fill") +  # position="fill"转换为比例
  labs(y = "Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right") +
  scale_fill_viridis_d(option = "D")  # 使用颜色渐变

# 标准化（按行或列）
heatmap_data_norm_4 <- prop.table(table(meta_df$L4_second_type, meta_df$layer), margin = 2)  # 按列标准化

# 绘制热图
p<-pheatmap(heatmap_data_norm_4,
            color = colorRampPalette(c("white", "red"))(100),
            cluster_rows = TRUE,
            cluster_cols = TRUE,
            #main = "Frequency of L4_second_type by mp",
            angle_col = '45',
            fontsize_row = 15,
            fontsize_col = 15)
ggsave('Frequency of L4_second_type by layer.pdf',plot = p,width = 6, height = 12, units = "in")
# 计算交叉频数表L3
freq_table_3 <- as.data.frame(table(meta_df$L3_second_type, meta_df$layer))
colnames(freq_table_3) <- c("L3_type", "layer", "Count")

# 绘制堆叠条形图
ggplot(freq_table_3, aes(x =L3_type , y = Count, fill = layer)) +
  geom_bar(stat = "identity", position = "fill") +  # position="fill"转换为比例
  labs(y = "Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right") +
  scale_fill_viridis_d(option = "D")  # 使用颜色渐变

# 标准化（按行或列）
heatmap_data_norm_3 <- prop.table(table(meta_df$L3_second_type, meta_df$layer), margin = 2)  # 按列标准化

# 绘制热图
pheatmap(heatmap_data_norm_3,
         color = colorRampPalette(c("white", "red"))(100),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         #main = "Frequency of L4_second_type by mp",
         angle_col = '45',
         fontsize_row = 15,
         fontsize_col = 15)




# 计算交叉频数表L2
freq_table_2 <- as.data.frame(table(meta_df$L2_second_type, meta_df$layer))
colnames(freq_table_2) <- c("L2_type", "layer", "Count")

# 绘制堆叠条形图
ggplot(freq_table_2, aes(x =L2_type , y = Count, fill = layer)) +
  geom_bar(stat = "identity", position = "fill") +  # position="fill"转换为比例
  labs(y = "Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right") +
  scale_fill_viridis_d(option = "D")  # 使用颜色渐变

# 标准化（按行或列）
heatmap_data_norm_2 <- prop.table(table(meta_df$L2_second_type, meta_df$layer), margin = 2)  # 按列标准化

# 绘制热图
pheatmap(heatmap_data_norm_2,
         color = colorRampPalette(c("white", "red"))(100),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         #main = "Frequency of L4_second_type by mp",
         angle_col = '45',
         fontsize_row = 15,
         fontsize_col = 15)


