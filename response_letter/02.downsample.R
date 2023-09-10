library(Seurat)
library(ggsci)
library(ggrepel)
library(dplyr)
col <- c(pal_npg("nrc")(10)[1:6],"#FAFD7CFF","#FF6F00FF",pal_lancet("lanonc")(9)[c(1,3,7)],"#660099FF","#B5CF6BFF","#B24745FF","#CCFF00FF",
         "#FFCD00FF","#800000FF","#20854EFF","#616530FF","#FF410DFF","#EE4C97FF","#FF1463FF","#00FF00FF","#990080FF","#00FFFFFF",
         "#666666FF","#CC33FFFF","#00D68FFF","#4775FFFF","#C5B0D5FF","#FDAE6BFF","#79CC3DFF","#996600FF","#FFCCCCFF","#0000CCFF",
         "#7A65A5FF","#1A5354FF","#24325FFF")

#ccle <- readRDS("/Users/shangliu/01.terms/01.lung_metastasis/05.Other/new/cellline5_20210727.rds")
ccle <- readRDS("/Users/shangliu/01.terms/03.ecDNA/NC/01.Figure1/Hep1-6_cellline.integrated.integrated.Find.rds")
DimPlot(ccle)
ann <- ccle@meta.data

table(ann$Batch)

pos <- which(ccle$Batch %in% "NTC_cellline")
cellname <- colnames(ccle)[pos]
mat <- ccle@assays$RNA@counts[,pos]

pbmc <- CreateSeuratObject(mat)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc,pattern = "^mt")

pbmc <- subset(pbmc,percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features =   VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.2)
pbmc <- RunUMAP(pbmc, dims = 1:10)

#VlnPlot(pbmc,features = "percent.mt")

DimPlot(pbmc,label = F,cols = col)
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc, features = top10$gene) 
pos <- which(pbmc.markers$cluster %in% 2)
pbmc.markers1 <- pbmc.markers[pos,]
#FeaturePlot(pbmc,features = "percent.mt")

#############
set.seed(2)
cellname <- sample(colnames(pbmc),300)
pbmc.sub <- subset(pbmc,cells = cellname)
pbmc.sub$old.clusters <- pbmc.sub$seurat_clusters

pbmc.sub <- NormalizeData(pbmc.sub, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc.sub <- FindVariableFeatures(pbmc.sub, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc.sub)
pbmc.sub <- ScaleData(pbmc.sub, features = all.genes)
pbmc.sub <- RunPCA(pbmc.sub, features =   VariableFeatures(object = pbmc.sub))
pbmc.sub <- FindNeighbors(pbmc.sub, dims = 1:10)
pbmc.sub <- FindClusters(pbmc.sub, resolution = 0.45)
pbmc.sub <- RunUMAP(pbmc.sub, dims = 1:10)

DimPlot(pbmc.sub,label = T,cols = col,group.by = "old.clusters")
DimPlot(pbmc.sub,label = T,cols = col)

table(pbmc.sub$old.clusters,pbmc.sub$seurat_clusters)

#pos <- which(pbmc.sub$seurat_clusters %in% 2)
#pbmc.sub$seurat_clusters[pos] <- 0

Idents(pbmc.sub) <- factor(pbmc.sub$seurat_clusters,levels = c(0,1,2,3,4))
DoHeatmap(pbmc.sub, features = top10$gene) 





