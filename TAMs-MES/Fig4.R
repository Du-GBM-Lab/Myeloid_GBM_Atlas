knitr::opts_chunk$set(echo = TRUE)


## Figure 2D-G - Co-localisation analysis of macrophage subsets, glioblastoma cellular subtype and anatamical niches

library(Seurat)
library(nichenetr)
library(readxl)
library(dplyr)
library(readr)
library(nichenetr)
library(dplyr)
library(corrplot)
library(PerformanceAnalytics)
library(RColorBrewer)
library(ggplot2)
# Load data
gbm.merged <- readRDS("D:/ProgramData/Spatial Gene Expression/Macrophage-mediated myelin recycling fuels brain cancer malignancy/GEO Data processed/visium/visium_merged_complete.rds")
colnames(gbm.merged@meta.data)
DimPlot(gbm.merged,group.by = 'Cell.Subtype')
DimPlot(gbm.merged,group.by = 'celltype')
DimPlot(gbm.merged,group.by = 'TME')
DimPlot(gbm.merged,group.by = 'celltype')
DimPlot(gbm.merged,group.by = 'MDM_subset')
DimPlot(gbm.merged,group.by = 'MG_subset')
DimPlot(gbm.merged,group.by = 'orig.ident')
DimPlot(gbm.merged,group.by = 'Location')
DimPlot(gbm.merged,group.by = 'GPNMBhigh_sig')
DimPlot(gbm.merged,group.by = 'Location')
DimPlot(gbm.merged,group.by = 'Location')
DimPlot(gbm.merged,group.by = 'Location')

table(gbm.merged$Cell.Subtype,gbm.merged$celltype)


## Figure 2E-G  P53_Rec_S4_2
P1<-SpatialDimPlot(gbm.merged[, gbm.merged$TME %in% c("MG1-P2RY12",	"MG2-TNF", "MG3-GPNMB" ,	"MDM1-CCR7",	"MDM2-H2-EB1", "MDM3-GPNMB")], 
                   group.by = "GPNMBhigh_sig",   stroke = 0,   image.alpha = 1, alpha = 1, 
                   images = c("Ink4a_Prim_S3")) + scale_fill_manual(values =  c("#FFA07A", "#AB82FF"))
P2<-SpatialDimPlot(gbm.merged, group.by = "Cell.Subtype",stroke = 0,  images = c("Ink4a_Prim_S3"), 
                   image.alpha = 0) + scale_fill_manual(values = c("darkolivegreen3", "pink", "tomato3","darkgoldenrod1", "darkgoldenrod3", "royalblue3", "NA" ))
P3<-SpatialDimPlot(gbm.merged, group.by = "Location",stroke = 0,  images = c("Ink4a_Prim_S3"), 
                   image.alpha =  0) + scale_fill_manual(values = c("orange", "violetred1", "darkolivegreen3","lightblue", "grey")) 
P1+P2+P3
## Supplementary Figure 7C-E (additional representitive images of primary tumor)
P1<-SpatialDimPlot(gbm.merged[, gbm.merged$TME %in% c("MG1-P2RY12",	"MG2-TNF", "MG3-GPNMB" ,	"MDM1-CCR7",	"MDM2-H2-EB1", "MDM3-GPNMB")], 
                   group.by = "GPNMBhigh_sig",   stroke = 0,   image.alpha = 1, alpha = 1, 
                   images = c("P53_Prim_S1_2")) + scale_fill_manual(values =  c("#FFA07A", "#AB82FF"))
P2<-SpatialDimPlot(gbm.merged, group.by = "Cell.Subtype",stroke = 0,  images = c("P53_Prim_S1_2"), 
                   image.alpha = 0) + scale_fill_manual(values = c("darkolivegreen3", "pink", "tomato3","darkgoldenrod1", "darkgoldenrod3", "royalblue3", "NA" ))
P3<-SpatialDimPlot(gbm.merged, group.by = "Location",stroke = 0,  
                   images = c("P53_Prim_S1_2"), image.alpha =  0) + scale_fill_manual(values = c("orange", "violetred1", "darkolivegreen3","lightblue", "grey")) 
P1+P2+P3

SpatialDimPlot(gbm.merged,ncol = 3)

# Figure 2E-G  P53_Rec_S4_2
P1<-SpatialDimPlot(gbm.merged[, gbm.merged$TME %in% c("MG1-P2RY12",	"MG2-TNF", "MG3-GPNMB" ,	"MDM1-CCR7",	"MDM2-H2-EB1", "MDM3-GPNMB")], 
                   group.by = "GPNMBhigh_sig",   stroke = 0,   image.alpha = 1, alpha = 1, 
                   images = c("P53_Rec_S4_2")) + scale_fill_manual(values =  c("#FFA07A", "#AB82FF"))
P2<-SpatialDimPlot(gbm.merged, group.by = "Cell.Subtype",stroke = 0,  images = c("P53_Rec_S4_2"), 
                   image.alpha = 0) + scale_fill_manual(values = c("darkolivegreen3", "tomato3", "tomato3","darkgoldenrod3", "darkgoldenrod3", "royalblue3", "NA" ))
P3<-SpatialDimPlot(gbm.merged, group.by = "Location",stroke = 0,  images = c("P53_Rec_S4_2"), 
                   image.alpha =  0) + scale_fill_manual(values = c("orange", "violetred1", "darkolivegreen3","lightblue", "grey")) 
P1+P2+P3

# Figure 2E-G Ink4a_Rec_S2
P1<-SpatialDimPlot(gbm.merged[, gbm.merged$TME %in% c("MG1-P2RY12",	"MG2-TNF", "MG3-GPNMB" ,	"MDM1-CCR7",	"MDM2-H2-EB1", "MDM3-GPNMB","NA")], 
                   group.by = "GPNMBhigh_sig",   stroke = 0,   image.alpha = 0, alpha = 1, 
                   images = c("Ink4a_Rec_S2")) + scale_fill_manual(values =  c("#FFA07A", "#AB82FF"))
P2<-SpatialDimPlot(gbm.merged, group.by = "Cell.Subtype",stroke = 0,  images = c("Ink4a_Rec_S2"), 
                   image.alpha = 0) + scale_fill_manual(values = c("darkolivegreen3", "tomato3", "tomato3","darkgoldenrod3", "darkgoldenrod3", "royalblue3", "NA" ))
P3<-SpatialDimPlot(gbm.merged, group.by = "Location",stroke = 0,  images = c("Ink4a_Rec_S2"), 
                   image.alpha =  0) + scale_fill_manual(values = c("royalblue3", "darkolivegreen3", "tomato3","darkgoldenrod3", "grey")) 
P1+P2+P3

table(gbm.merged$TME)

c( "darkolivegreen3", "tomato3", "tomato3","darkgoldenrod3", "darkgoldenrod3", "royalblue3" ))























gbm.merged_tumor <- gbm.merged[, gbm.merged$tissue %in% c("Tumor")]
gbm.merged_tumormac <- gbm.merged_tumor[, gbm.merged_tumor$Macrophage_score %in% c("Macrophage")]

Idents(gbm.merged_tumormac) <- gbm.merged_tumormac$orig.ident
DimPlot(gbm.merged_tumormac)
gbm.merged_tumormac <- RenameIdents(object = gbm.merged_tumormac, "S1_2" = "S9",
                                    "S2_2" = "S10",
                                    "S3_2" = "S11",
                                    "S4_2" = "S12",
                                    "S8_2" = "S16"
)
gbm.merged_tumormac$orig.ident <- Idents(gbm.merged_tumormac)
correlation <- gbm.merged_tumormac@meta.data


## Make LLM score
DefaultAssay(gbm.merged_tumormac) <- "SCT"

## GPNMBhigh gene signature
Gene <- c("Gpnmb",  "Fabp5",  "Hmox1",  "Spp1" ,  "Arg1")

##Add score to every cell and create new seurat object. Not necessary to create a new object but in case anything goes wrong your original data is still there.
##This part is taken from the cellcycle classification, however I only used the part where Geneset scores were assigned to the cell
##create a score for every group

gbm.merged_tumormac <- NormalizeData(gbm.merged_tumormac) # Needs to be normalised 
fdscores <- AddModuleScore(gbm.merged_tumormac, features= list(c(Gene)), name="Gene",nbin=100)

##define the different groups in your Genelist, add '1' to every groupname. Groupnames can be checked in the metadata -> str(fdscores@meta.data)
groups <- c("Gene1")

##load function to create density values
densMode <- function(x){
  td <- density(x)
  tdx <- td$x
  tdy <- td$y
  minx <- tdx[which(diff(sign(diff(tdy)))==2)]
  peakx <- tdx[which(diff(sign(diff(tdy)))==-2)]
  return(list(minx=minx, maxy=peakx))
}

##For every group determine the thesholds and plot several plots
for (i in groups){
  ##create densityplots
  plot(density(fdscores@meta.data[,i]), main=paste("densityplot"))
  ##classify the cells based on thresholds of 0.1
  gbm.merged@meta.data[,paste("assignedto",i, sep="")] <- "nonclassified"
  gbm.merged@meta.data[which(fdscores@meta.data[,i]>0.1),paste("assignedto",i, sep="")] <- paste("assignedto",i, sep=" ")
  fdscores_gpnmbhigh <- fdscores
}

correlation_gpnmbhigh <- fdscores_gpnmbhigh@meta.data # acquired from Figure 1k

correlation_perwell <- data.frame(correlation$MES1.Score, correlation$MES2.Score, correlation$AC.Score, correlation$NPC1.Score, correlation$NPC2.Score, correlation$OPC.Score,  correlation$PAN.Score, correlation$LE.Score, correlation$MVP.Score, correlation$CT.Score,  correlation$M01.Score, correlation$M02.Score, correlation$M03.Score,correlation$B01.Score, correlation$B02.Score, correlation$B03.Score, correlation_gpnmbhigh$Gene1)


corrplot(cor(correlation_perwell[-c(631),]),  
         tl.col = "black", # Labels color
         bg = "white",     # Background color
         col=brewer.pal(n=8, name="PRGn"))       # Color palette

chart.Correlation(cor(correlation_perwell[-c(631),]), histogram=TRUE, pch=20)

write.csv2(cor(correlation_perwell[-c(631),]), "~/Desktop/Correlation_Figure2D.csv")
rm(gbm.merged_tumor)
rm(gbm.merged_tumormac)



# Figure 2D (2): Assessing the fold-change of the various cell types and tumor areas between primary and recurrent PDG-Ink4a and PDG-shp53 mice.

gbm.merged_quantify <- gbm.merged[, gbm.merged$tissue %in% c("Tumor")]
Idents(gbm.merged_quantify) <- gbm.merged_quantify$orig.ident
gbm.merged_quantify <- RenameIdents(object = gbm.merged_quantify, "S1_2" = "S9",
                                    "S2_2" = "S10",
                                    "S3_2" = "S11",
                                    "S4_2" = "S12",
                                    "S8_2" = "S16"
)
gbm.merged_quantify$orig.ident <- Idents(gbm.merged_quantify)

table_TME <- table(gbm.merged_quantify$TME, gbm.merged_quantify$orig.ident)
table_CellSubtype <- table(gbm.merged_quantify$Cell.Subtype, gbm.merged_quantify$orig.ident)
table_IVY <- table(gbm.merged_quantify$Location, gbm.merged_quantify$orig.ident)
table_GPNMBhigh <- table(gbm.merged_quantify$GPNMBhigh_sig, gbm.merged_quantify$orig.ident)

library(openxlsx)
write.xlsx(table_TME, file="~/Desktop/table_fractions-TME.xlsx", sheetName="table_TME", append=TRUE, rowNames = F)
write.xlsx(table_CellSubtype, file="~/Desktop/table_fractions-Neftel.xlsx", sheetName="table_CellSubtype", append=TRUE, rowNames = F)
write.xlsx(table_IVY, file="~/Desktop/table_fractions-IVY.xlsx", sheetName="table_IVY", append=TRUE, rowNames = F)
write.xlsx(table_GPNMBhigh, file="~/Desktop/table_fractions-GPNMBhigh.xlsx", sheetName="table_GPNMBhigh", append=TRUE, rowNames = F)



# Figure 2E-G Representative visualization of VISIUM 10X spatial transcriptomic analyses in recurrent PDG-Ink4a glioblastoma, highlighting GPNMBhigh deserted and enriched areas

```{r}
Idents(gbm.merged) <- gbm.merged$TME
my_levels <-  c("MG1-P2RY12",	"MG2-TNF", "MG3-GPNMB" ,	"MDM1-CCR7",	"MDM2-H2-EB1", "MDM3-GPNMB", "NA")

# Relevel object@ident
Idents(gbm.merged) <- factor(Idents(gbm.merged), levels = my_levels)
gbm.merged$TME <- Idents(gbm.merged)

ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

color_list <- ggplotColours(n=6)
color_list <- c(color_list, "grey")

Idents(gbm.merged) <- gbm.merged$TME

SpatialDimPlot(gbm.merged, group.by = "TME",   stroke = 0,   image.alpha = 0, alpha = 1, 
               images = c(" P53_Prim_S1_2")) + scale_fill_manual(values = color_list)

SpatialDimPlot(gbm.merged,  stroke = 0,   image.alpha = 0, alpha = 1, images = c("P53_Prim_S1_2"), cells.highlight = CellsByIdentities(gbm.merged), facet.highlight = TRUE,  ncol = 4)
table(gbm.merged$TME)
table(gbm.merged$celltype)
table(gbm.merged$Location)
table(gbm.merged$Cell.Subtype)
table(gbm.merged$TME,gbm.merged$Cell.Subtype)
names(gbm.merged)
