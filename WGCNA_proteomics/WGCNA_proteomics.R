library(tidyverse)
library(jsonlite)
library(WGCNA)
library(flashClust)
library(matrixStats)

enableWGCNAThreads()
# Set number of threads
nThreads = 15

# Set directlory for storing experiment results
setwd("d:/ProgramData/Glioma-immunelncRNA/WGCNA_proteomics/")

getwd()
vst_mat = read.csv("processed/lfq.csv", row.names=1)
anno = read.csv("./processed/metadata.csv", row.names=1)
vst_mat_corrected <- empiricalBayesLM(t(vst_mat), anno$batch)
vst_mat <- t(vst_mat_corrected[[1]])
anno$group <- anno[, "Mey.subclass..0.low..1.high."]
conditions<-factor(sapply(rownames(anno),function(id){paste(anno[id,c("group")],collapse = '_')}))
conditions
trait_df<-data.frame(t(sapply(conditions,function(condition){table(condition)})),row.names=rownames(anno))
# This creates an object called "datExpr" that contains the normalized counts file output from DESeq2
datExpr = vst_mat
datExpr = t(datExpr)
colnames(trait_df)
colnames(trait_df) <- c("low", "high")
datTraits = trait_df
table(rownames(datTraits)==rownames(datExpr))



# Cluster samples by expression ----------------------------------------------------------------

A = adjacency(t(datExpr),type="signed") # this calculates the whole network connectivity
k = as.numeric(apply(A,2,sum))-1 # standardized connectivity
Z.k = scale(k)
thresholdZ.k = -2.5 # often -2.5
outlierColor = ifelse(Z.k<thresholdZ.k,"red","black")
sampleTree = flashClust(as.dist(1-A), method = "average")
# Convert traits to a color representation where red indicates high values
traitColors = data.frame(numbers2colors((datTraits),signed=FALSE))
dimnames(traitColors)[[2]] = paste(names(datTraits))
datColors = data.frame(outlier = outlierColor,traitColors)

pdf(filename="sampleDendorgramAndTraitHeatmap.pdf")
plotDendroAndColors(sampleTree,groupLabels=names(datColors),
                    colors=datColors,main="Sample Dendrogram and Trait Heatmap")
dev.off()

##remove outliers
outliers <- c("F17","F1")
gsg = goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK
##Flagging genes and samples with too many missing values.....step 1
datTraits <- datTraits[!(row.names(datTraits) %in% outliers),]
datExpr <- datExpr[!(row.names(datExpr) %in% outliers),]
table(rownames(datTraits)==rownames(datExpr))
# Choose a soft threshold power
powers = c(c(1:10), seq(from =10, to=100, by=0.5)) #choosing a set of soft-thresholding powers
sft = pickSoftThreshold(datExpr, dataIsExpr = TRUE, powerVector=powers, 
                        verbose =5, networkType="signed",
                        corOptions = list(use = 'p')) #call network topology analysis function

sizeGrWindow(9,5)
par(mfrow= c(1,2))
cex1=0.9
png("scaleIndependence.png")
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab= "Soft Threshold (power)", 
     ylab="Scale Free Topology Model Fit, signed R^2", type= "n", main= paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=cex1, col="red")
abline(h=0.6, col="red")
dev.off()

png("meanConnectivity.png")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab= "Soft Threshold (power)", ylab="Mean Connectivity", type="n", 
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")
dev.off()
softPower <- 9
nThreads <- 15
enableWGCNAThreads(nThreads=nThreads)
adjacency = adjacency(datExpr, power = softPower, type = "signed") #specify network type
# Construct Networks- USE A SUPERCOMPUTER IRL -----------------------------

#translate the adjacency into topological overlap matrix and calculate the corresponding dissimilarity:
TOM = TOMsimilarity(adjacency, TOMType="signed") # specify network type
dissTOM = 1-TOM

# Generate Modules --------------------------------------------------------

# Generate a clustered gene tree
geneTree = flashClust(as.dist(dissTOM), method="average")

png("geneClusteringTOM.png")
plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE, hang=0.04)
dev.off()
packageVersion("clusterProfiler")
#This sets the minimum number of genes to cluster into a module
minModuleSize = 10
x = 4 
dynamicMods = cutreeDynamic(dendro = geneTree, distM = as.matrix(dissTOM), 
                            method="hybrid", pamStage = F, deepSplit = x, 
                            minClusterSize = minModuleSize)

dynamicColors= labels2colors(dynamicMods)
MEList= moduleEigengenes(datExpr, colors= dynamicColors)#,softPower = 14)
MEs= MEList$eigengenes
MEDiss= 1-cor(MEs)
METree= flashClust(as.dist(MEDiss), method= "average")

save(dynamicMods, MEList, MEs, MEDiss, METree, file= "Network_allSamples_signed_RLDfiltered.RData")

#plots tree showing how the eigengenes cluster together
png(paste("eigenGenesClustering", ".png", sep=""))
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")
dev.off()

#softPower = 10 # see scale independence file "softthreshold_SOD1.png"

#### -----------------------------------------------------------------------
#set a threhold for merging modules. In this example we are not merging so MEDissThres=0.0
MEDissThres = 0.4
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight= MEDissThres, verbose =3)
mergedColors = merge$colors
mergedMEs = merge$newMEs

#plot dendrogram with module colors below it
png("geneTree.png")
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)
dev.off()
moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

save(MEs, moduleLabels, moduleColors, geneTree, file= "Network_allSamples_signed_nomerge_RLDfiltered.RData")

#Define number of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use="p")

write.csv(moduleTraitCor, "moduleTraitCor.csv")

moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

write.csv(moduleTraitPvalue, "moduleTraitPvalue.csv")
#display the corelation values with a heatmap plot 我自己改了一下下下
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)


textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar=  c(6, 9, 3, 3))
png("moduleTraitRelationshipHeatmap_corrAll.png", width=600, height=600)
labeledHeatmap(Matrix = t(moduleTraitCor),
               yLabels = names(datTraits),
               xLabels = names(MEs),
               xSymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = t(textMatrix),
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

#get selected genes
selected_Colors<-sapply(names(MEs),function(ME){sub("ME","",ME)})
#names(datExpr)[moduleColors=="brown"]
module_genes<-sapply(selected_Colors,function(module_color){rownames(t(datExpr)[mergedColors==module_color,])})

WGCNA_folder<-"./"
saveRDS(module_genes,file=file.path(WGCNA_folder,'module_genes.rds'))
write(toJSON(module_genes),file=file.path(WGCNA_folder,'module_genes.json'))

## enrichment of each module
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(DOSE)
  library(org.Mm.eg.db)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(pathview)
  library(tidyverse)
})

make_go_enr = function(genelist, 
                       fname, 
                       folder, 
                       fromType="SYMBOL",
                       #OrgDb=org.Mm.eg.db,
                       OrgDb=org.Hs.eg.db,
                       prefixes=c("BP","MF","CC","All"),
                       image_width=30,
                       image_height=30){
  
  gene.df <- bitr((genelist), fromType = fromType,
                  toType = c("ENTREZID"),
                  OrgDb = OrgDb)
  for (prefix in prefixes) {
    ego <- enrichGO(gene          = gene.df$ENTREZID,
                    OrgDb         = OrgDb,
                    ont           = prefix,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05,
                    readable = T)
    
    dir.create(file.path(folder, fname))
    
    if (!is.null(ego) && nrow(ego)!=0){
      write.csv(ego, file.path(folder, fname, paste0("go_enr_", prefix,"_",fname,".csv")))
      
      p1 = dotplot(ego, showCategory=20) + ggtitle(paste0("GO enrichment ", prefix))
      ggsave(p1, filename = file.path(folder,fname, paste0("dotplot_go_", prefix,"_",fname,".pdf")), width = image_width, 
             height=image_height, units="cm")
      
      p2 = barplot(ego, showCategory=20) + ggtitle(paste0("GO enrichment ", prefix))
      ggsave(p2, filename = file.path(folder, fname, paste0("barplot_go_", prefix,"_",fname,".pdf")), width = image_width, 
             height=image_height, units="cm")
    }
  }
}
make_go_enr1 = function(genelist, 
                        fname, 
                        folder, 
                        fromType="SYMBOL",
                        #OrgDb=org.Mm.eg.db,
                        OrgDb=org.Hs.eg.db,
                        prefixes=c("BP","MF","CC","All"),
                        image_width=20,
                        image_height=20){
  
  gene.df <- bitr((genelist), fromType = fromType,
                  toType = c("ENTREZID"),
                  OrgDb = OrgDb)
  for (prefix in prefixes) {
    ego <- enrichGO(gene          = gene.df$ENTREZID,
                    OrgDb         = OrgDb,
                    ont           = prefix,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05,
                    readable = T)
    
    dir.create(file.path(folder, fname))
    
    if (!is.null(ego) && nrow(ego)!=0){
      write.csv(ego, file.path(folder, fname, paste0("go_enr_", prefix,"_",fname,".csv")))
      
      p1 = dotplot(ego, showCategory=10, font.size=20) + ggtitle(paste0("GO enrichment ", prefix))
      ggsave(p1, filename = file.path(folder,fname, paste0("dotplot_go_", prefix,"_",fname,".pdf")), width = image_width, 
             height=image_height, units="cm")
      
      p2 = barplot(ego, showCategory=10, font.size=20) + ggtitle(paste0("GO enrichment ", prefix))
      ggsave(p2, filename = file.path(folder, fname, paste0("barplot_go_", prefix,"_",fname,".pdf")), width = image_width, 
             height=image_height, units="cm")
    }
  }
}
modules_ <- c("MEbrown", "MEdarkorange")
result_folder<-"results_MEbrown_MEdarkorange"
overview_df<- data.frame()

dir.create(result_folder)

for (module in modules_){ 
  
  make_go_enr1(module_genes[[module]],module,file.path(WGCNA_folder,result_folder))
}

result_folder<-"results"
overview_df<- data.frame()

dir.create(result_folder)

for (module in names(module_genes)){ 
  
  make_go_enr(module_genes[[module]],module,file.path(WGCNA_folder,result_folder))
}








##