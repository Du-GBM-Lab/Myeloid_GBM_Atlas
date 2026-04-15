library(ggplot2)
library(dplyr)
library(ggrepel)
library(RColorBrewer)
library(Seurat)

Core <- qread("D:/ProgramData/GBMAP-data/Core.GBmap-u-t.qs") 
colnames(Core@meta.data)
# 提取author的分布
author_counts <- as.data.frame(table(Core$author))
colnames(author_counts) <- c("Author", "Count")

# 计算百分比
author_counts$Percentage <- round(author_counts$Count / sum(author_counts$Count) * 100, 1)

# 按数量排序
author_counts <- author_counts[order(-author_counts$Count), ]
# 获取作者类别的数量
author_count <- length(unique(Core$author))

# 创建扩展调色板
# 从Set3调色板的基础颜色开始，生成我们所需数量的颜色
getPalette <- colorRampPalette(brewer.pal(12, "Set3")) # Set3本身最多有12种颜色[8](@ref)
# 或者使用其他定性调色板作为基础，如Paired有12色，Set1有9色[6](@ref)

# 修改您的绘图代码，使用自定义的颜色标度
ggplot(author_counts, aes(x = "", y = Count, fill = Author)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme_void() +
  labs(title = paste("Total:", sum(author_counts$Count), "samples"),
       subtitle = paste("Number of authors:", nrow(author_counts))) +
  scale_fill_manual(values = getPalette(author_count)) + # 使用扩展后的颜色
  geom_text_repel(aes(label = paste0(Author, "\n", Count, " (", Percentage, "%)")),
                  position = position_stack(vjust = 0.5),
                  size = 3.5, max.overlaps = 20) +
  theme(plot.title = element_text(hjust = 0.8, face = "bold"),
        plot.subtitle = element_text(hjust = 0.8),
        legend.position = "none") # 如果类别太多，可以考虑隐藏图例


# 提取assay的分布
assay_counts <- as.data.frame(table(Core$assay))
colnames(assay_counts) <- c("Assay", "Count")

# 计算百分比
assay_counts$Percentage <- round(assay_counts$Count / sum(assay_counts$Count) * 100, 1)

# 按数量排序
assay_counts <- assay_counts[order(-assay_counts$Count), ]

# 获取assay类别的数量
assay_count <- length(unique(Core$assay))

# 创建扩展调色板
getPalette <- colorRampPalette(brewer.pal(12, "Set3"))

ggplot(assay_counts, aes(x = "", y = Count, fill = Assay)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme_void() +
  labs(title = paste("Total:", sum(assay_counts$Count), "samples"),
       subtitle = paste("Number of assays:", nrow(assay_counts))) +
  scale_fill_manual(values = getPalette(assay_count)) +
  geom_text_repel(aes(label = paste0(Assay, "\n", Count, " (", Percentage, "%)")),
                  position = position_stack(vjust = 0.5),
                  size = 3.5, max.overlaps = 20) +
  theme(plot.title = element_text(hjust = 0.8, face = "bold"),
        plot.subtitle = element_text(hjust = 0.8),
        legend.position = "none")


# 提取annotation_level_3的分布
assay_counts <- as.data.frame(table(Core$annotation_level_3))
colnames(assay_counts) <- c("Assay", "Count")

# 计算百分比
assay_counts$Percentage <- round(assay_counts$Count / sum(assay_counts$Count) * 100, 1)

# 按数量排序
assay_counts <- assay_counts[order(-assay_counts$Count), ]

# 获取assay类别的数量
assay_count <- length(unique(Core$annotation_level_3))

# 创建扩展调色板
getPalette <- colorRampPalette(brewer.pal(12, "Set3"))

ggplot(assay_counts, aes(x = "", y = Count, fill = Assay)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme_void() +
  labs(title = paste("Total:", sum(assay_counts$Count), "samples"),
       subtitle = paste("Number of assays:", nrow(assay_counts))) +
  scale_fill_manual(values = getPalette(assay_count)) +
  geom_text_repel(aes(label = paste0(Assay, "\n", Count, " (", Percentage, "%)")),
                  position = position_stack(vjust = 0.5),
                  size = 3.5, max.overlaps = 20) +
  theme(plot.title = element_text(hjust = 0.8, face = "bold"),
        plot.subtitle = element_text(hjust = 0.8),
        legend.position = "none")



