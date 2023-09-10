library(Seurat)
library(ggsci)
library(ggrepel)
library(dplyr)
library(ggplot2)

col <- c(pal_npg("nrc")(10)[1:6],"#FAFD7CFF","#FF6F00FF",pal_lancet("lanonc")(9)[c(1,3,7)],"#660099FF","#B5CF6BFF","#B24745FF","#CCFF00FF",
         "#FFCD00FF","#800000FF","#20854EFF","#616530FF","#FF410DFF","#EE4C97FF","#FF1463FF","#00FF00FF","#990080FF","#00FFFFFF",
         "#666666FF","#CC33FFFF","#00D68FFF","#4775FFFF","#C5B0D5FF","#FDAE6BFF","#79CC3DFF","#996600FF","#FFCCCCFF","#0000CCFF",
         "#7A65A5FF","#1A5354FF","#24325FFF")
###计算异质性分数
#2. 计算异质性分数

##compute diversity score
load("/Users/shangliu/01.terms/03.ecDNA/NC/01.Figure1/allcellines_update_newbatch_filter_anno.Rdata")

pca <- ob@reductions$pca@cell.embeddings[,1:30]
table(ob$Cell_line)
i <- "SKBR3"
temp_mat <- pca[which(ob$orig.ident==i),]
temp_mat1 <- data.frame(temp_mat)
ggplot(temp_mat1,aes(x=PC_1,y=PC_2))+geom_point() + theme_classic() 

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
temp_mat2 <- temp_mat1[colnames(pbmc),]
temp_mat2$clus <- pbmc$seurat_clusters
ggplot(temp_mat2,aes(x=PC_1,y=PC_2,color = clus))+geom_point() + theme_classic() 


###########
table(ob$Cell_line)
cellline <- "SNB75"
temp_mat <- pca[which(ob$orig.ident==cellline),]
temp_mat1 <- data.frame(temp_mat)
mat <- dist(temp_mat1)

mat <- as.matrix(mat)
siz <- NULL
for(i in 1:(nrow(mat)-1)){
  tmp <- mat[i,(i+1):ncol(mat)]
  siz <- c(siz,as.numeric(as.character(tmp)))
}
plot(density(siz),main = cellline)
pos <- which(siz > 5)
siz1 <- siz[-pos]
plot(density(siz1),main = cellline)


