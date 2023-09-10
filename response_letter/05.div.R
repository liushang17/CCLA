library(Seurat)
library(ggsci)
library(ggrepel)
library(dplyr)
col <- c(pal_npg("nrc")(10)[1:6],"#FAFD7CFF","#FF6F00FF",pal_lancet("lanonc")(9)[c(1,3,7)],"#660099FF","#B5CF6BFF","#B24745FF","#CCFF00FF",
         "#FFCD00FF","#800000FF","#20854EFF","#616530FF","#FF410DFF","#EE4C97FF","#FF1463FF","#00FF00FF","#990080FF","#00FFFFFF",
         "#666666FF","#CC33FFFF","#00D68FFF","#4775FFFF","#C5B0D5FF","#FDAE6BFF","#79CC3DFF","#996600FF","#FFCCCCFF","#0000CCFF",
         "#7A65A5FF","#1A5354FF","#24325FFF")
###计算异质性分数
#2. 计算异质性分数

##compute diversity score
load("/Users/shangliu/01.terms/03.ecDNA/NC/01.Figure1/allcellines_update_newbatch_filter_anno.Rdata")

pca <- ob@reductions$pca@feature.loadings
pca1 <- data.frame(pca)
pcgene <- NULL
for(i in 1:30){
  colnames(pca1)[i] <- "key"
  pca1 <- pca1[order(pca1$key),]
  genetmp1 <- rownames(pca1)[1:200]
  genetmp2 <- rownames(pca1)[1801:2000]
  tmp <- data.frame(pc = paste0("PC_",i), genename = c(genetmp1,genetmp2),type = rep(c("neg","pos"),each = 200))
  pcgene <- rbind(pcgene,tmp)
  colnames(pca1)[i] <- paste0("PC_",i)
}


i <- "SKBR3"
pos <- which(ob$orig.ident %in% i)
pbmc <- CreateSeuratObject(ob@assays$RNA@counts[,pos])
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc,pattern = "^mt")

#pbmc <- subset(pbmc,percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features =   VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:30)


DimPlot(pbmc,label = T,cols = col,reduction = "umap")  
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pos <- which(pbmc.markers$p_val_adj < 0.01)
pbmc.markers <- pbmc.markers[pos,]  
length(unique(pbmc.markers$gene))

geneall <- unique(as.character(pbmc.markers$gene))

pcinfo <- pbmc[['pca']]@feature.loadings

pos <- which(rownames(pcinfo) %in% rownames(pca1))

pbmcgene <- NULL
pca1 <- data.frame(pcinfo)
pbmcgene <- NULL
for(i in 1:30){
  colnames(pca1)[i] <- "key"
  pca1 <- pca1[order(pca1$key),]
  genetmp1 <- rownames(pca1)[1:50]
  genetmp2 <- rownames(pca1)[1951:2000]
  tmp <- data.frame(pc = paste0("PC_",i), genename = c(genetmp1,genetmp2),type = rep(c("neg","pos"),each = 50))
  pbmcgene <- rbind(pbmcgene,tmp)
  colnames(pca1)[i] <- paste0("PC_",i)
}

clus <- paste0("PC_",1:30)
siz <- NULL
for(i in 1:length(clus)){
  pos <- which(pcgene$pc %in% clus[i])
  pcgene1 <- pcgene[pos,]
  pos <- which(pcgene1$genename %in% rownames(pcinfo))
  siz <- c(siz,length(pos) / length(geneall))
}
summary(siz)
sort(siz)
siz



#mit <- data.frame(celline = c("Hela","SKBR3"),value = c(862,469))
mit <- data.frame(cellline = rep(c("Hela","SKBR3"),each = 30),value = c(tmpsiz,siz))
#ggplot(mit,aes(x=celline,y=value))+geom_bar(stat = "identity")+theme_classic()
ggplot(mit,aes(x=cellline,y = value))+geom_violin() +geom_boxplot() + theme_classic()
t.test(tmpsiz,siz,paired = T)


