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

pca <- ob@reductions$pca@cell.embeddings[,1:30]

i <- "Hela"
ITH <- c()
for (i in unique(ob$orig.ident)) {
  temp_mat <- pca[which(ob$orig.ident==i),]
  cutoff1 <- apply(temp_mat, 2, function(x) mean(x)-3*sd(x))
  cutoff2 <- apply(temp_mat, 2, function(x) mean(x)+3*sd(x))
  cell_bool <- apply(temp_mat, 1, function(x) sum(x[1:3] < cutoff1[1:3]) < 3 & sum(x[1:3] > cutoff2[1:3]) < 3)
  temp_mat <- temp_mat[cell_bool,]
  mean_m<-apply(temp_mat, 2, function(y) mean(y))
  diversity <- mean(apply(temp_mat, 1, function(x) sqrt(sum((x-mean_m)^2))))
  diversity
  ITH <- c(ITH,diversity)
  
}
info <- data.frame(cellline = unique(ob$orig.ident),value = ITH)

pca <- ob@reductions$umap@cell.embeddings

i <- "Hela"
ITH <- c()
for (i in unique(ob$orig.ident)) {
  temp_mat <- pca[which(ob$orig.ident==i),]
  cutoff1 <- apply(temp_mat, 2, function(x) mean(x)-3*sd(x))
  cutoff2 <- apply(temp_mat, 2, function(x) mean(x)+3*sd(x))
  cell_bool <- apply(temp_mat, 1, function(x) sum(x[1:2] < cutoff1[1:2]) < 3 & sum(x[1:2] > cutoff2[1:2]) < 3)
  temp_mat <- temp_mat[cell_bool,]
  mean_m<-apply(temp_mat, 2, function(y) mean(y))
  diversity <- mean(apply(temp_mat, 1, function(x) sqrt(sum((x-mean_m)^2))))
  diversity
  ITH <- c(ITH,diversity)
  
}
info$umap <- ITH
cor.test(info$umap,info$value)
ggplot(info,aes(x=value,y=umap))+geom_point()+geom_smooth(method = "lm") + theme_classic()

library(ggplot2)
i <- "SNB75"
temp_mat <- pca[which(ob$orig.ident==i),]
temp_mat1 <- data.frame(temp_mat)
ggplot(temp_mat1,aes(x=PC_1,y=PC_2))+geom_point() + theme_classic() 

pca <- ob[['pca']]@feature.loadings





denum <- NULL 
#ITH <- c()
for (i in unique(ob$orig.ident)) {
#i <- "SKBR3"
pos <- which(ob$orig.ident %in% i)
pbmc <- CreateSeuratObject(ob@assays$RNA@counts[,pos])
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc,pattern = "^mt")

#pbmc <- subset(pbmc,percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features =   VariableFeatures(object = pbmc))
#pbmc <- FindNeighbors(pbmc, dims = 1:30)
#pbmc <- FindClusters(pbmc, resolution = 0.5)
#pbmc <- RunUMAP(pbmc, dims = 1:30)

temp_mat <- pbmc@reductions$pca@cell.embeddings[,1:30]
cutoff1 <- apply(temp_mat, 2, function(x) mean(x)-3*sd(x))
cutoff2 <- apply(temp_mat, 2, function(x) mean(x)+3*sd(x))
cell_bool <- apply(temp_mat, 1, function(x) sum(x[1:3] < cutoff1[1:3]) < 3 & sum(x[1:3] > cutoff2[1:3]) < 3)
temp_mat <- temp_mat[cell_bool,]
mean_m<-apply(temp_mat, 2, function(y) mean(y))
diversity <- mean(apply(temp_mat, 1, function(x) sqrt(sum((x-mean_m)^2))))

denum <- c(denum,diversity)
#DimPlot(pbmc,label = T,cols = col,reduction = "umap")  
#pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#pos <- which(pbmc.markers$p_val_adj < 0.01)
#pbmc.markers <- pbmc.markers[pos,]  
#denum <- c(denum,length(unique(pbmc.markers$gene)))
}
mit <- data.frame(ith = ITH,denum = denum)
rownames(mit) <- unique(ob$orig.ident)
ggplot(mit,aes(x=ith,y = denum))+geom_point()+ geom_smooth(method = "lm")
cor.test(mit$ith,mit$denum)

geneall <- unique(as.character(pbmc.markers$gene))

pcinfo <- pbmc[['pca']]@feature.loadings

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


siz <- NULL
for(i in 1:length(clus)){
  pos <- which(pcgene$pc %in% clus[i])
  pcgene1 <- pcgene[pos,]
  pos <- which(pcgene1$genename %in% geneall)
  siz <- c(siz,length(pos) / length(geneall))
}
summary(siz)
sort(siz)
siz

pos <- which(geneall %in% pcgene$genename)
length(geneall)

clus <- paste0("PC_",1:30)
sim.scor <- NULL
for(i in 1:length(clus)){
  pos <- which(pbmcgene$pc %in% clus[i])
  pbmcgene1 <- pbmcgene[pos,]
  siz <- NULL
  for(j in 1:length(clus)){
    pos <- which(pcgene$pc %in% clus[j])
    pcgene1 <- pcgene[pos,]
    pos <- which(pcgene1$genename %in% pbmcgene1$genename)
    tmp <- length(pos) / length(unique(c(as.character(pcgene1$genename),as.character(pbmcgene1$genename))))
    siz <- c(siz,tmp)
  }
  tmp1 <- data.frame(pc = clus[i],score = siz)
  sim.scor <- c(sim.scor,max(siz))
}

summary(sim.scor)

  #denum <- c(denum,length(unique(pbmc.markers$gene)))
#}  
#info$denum <- denum
info$new_value <- ITH
  ggplot(info,aes(x=new_value,y=log2(denum)))+geom_point()+ labs(x="diveristy score", y = "DEG number")  
cor.test(info$new_value,(info$denum))  
  
pcinfo <- Loadings(ob[["pca"]])
saveRDS(pcinfo,file = "/Users/shangliu/01.terms/03.ecDNA/NC/01.Figure1/pcinfo.rds")



  