library(Seurat)
library(dplyr)
t1 <- readRDS("/Users/shangliu/01.terms/03.ecDNA/NC/07.hy/T1_MDAMB231.rds")
t2 <- readRDS("/Users/shangliu/01.terms/03.ecDNA/NC/07.hy/T2_MDAMB231.rds")

t1.ann <- read.table("/Users/shangliu/01.terms/03.ecDNA/NC/07.hy/MDAMB231_T1.fixed_SW.txt",sep = "\t",header = F)
t2.ann <- read.table("/Users/shangliu/01.terms/03.ecDNA/NC/07.hy/MDAMB231_T2.fixed_SW.txt",sep = "\t",header = F)
t1.ann <- unique(t1.ann)
t2.ann <- unique(t2.ann)
#############
pos <- which(t1.ann$V1 %in% 0)
t1.ann1 <- t1.ann[pos,]
t1.ann2 <- t1.ann[-pos,]
pos <- which(t1.ann1$V2 %in% t1.ann2$V2)
t1.ann3 <- t1.ann1[-pos,]
sui <- data.frame(table(t1.ann3$V2))

pos <- which(t1.ann1$V2 %in% "GACTGTGTCTCTGACAGTCTCAGAGTGT")
t1.ann4 <- t1.ann1[pos,]
pos <- which(t2.ann$V2 %in% "GACTGTGTCTCTGACAGTCTCAGAGTGT") #CTGTGACACTCAGTGTGTCTGTCACACA GACTGTGACTCTGACAGTGAGTCTGTCT
t2.ann2 <- t2.ann[pos,]
table(t2.ann2$V1)

pos <- which(colnames(t1) %in% t1.ann4$V3)
mat <- data.frame(as.matrix(t1@assays$RNA@data[,pos]))
pos <- (which(colnames(t2) %in% t2.ann2$V3))
mit <- data.frame(as.matrix(t2@assays$RNA@data[,pos]))

t2.ann2 <- t2.ann2[order(t2.ann2$V1),]
mit <- mit[,as.character(t2.ann2$V3)]
pos <- which(rownames(mat) %in% rownames(mit))
mat1 <- mat[pos,]
mit1 <- mit[rownames(mat1),]
mat2 <- cbind(mat1,mit1)
ann <- data.frame(clus = c(rep("T1",ncol(mat1)),as.character(t2.ann2$V1)))
rownames(ann) <- colnames(mat2)
#pheatmap::pheatmap(cor(mat2),cluster_rows = F,cluster_cols = F,annotation_row = ann,annotation_col = ann)

clus <- unique(ann$clus)
mitt <- matrix(nrow = nrow(mat2),ncol = length(clus))
for(i in 1:length(clus)){
  pos <- which(ann$clus %in% clus[i])
  matt <- mat2[,pos]
  mitt[,i] <- rowMeans(matt)
}
colnames(mitt) <- clus
pheatmap::pheatmap(cor(mitt),show_rownames = T,show_colnames = F,cluster_rows = F,cluster_cols = F)


###############
pos <- which(t1.ann1$V2 %in% "GTGTGTCTGAGAGAGTCTCTGACAGAGA")
t1.ann4 <- t1.ann1[pos,]
pos <- which(t2.ann$V2 %in% "GTGTGTCTGAGAGAGTCTCTGACAGAGA") #CTGTGACACTCAGTGTGTCTGTCACACA GACTGTGACTCTGACAGTGAGTCTGTCT
t2.ann2 <- t2.ann[pos,]
table(t2.ann2$V1)

pos <- which(colnames(t1) %in% t1.ann4$V3)
mat <- data.frame(as.matrix(t1@assays$RNA@data[,pos]))
pos <- (which(colnames(t2) %in% t2.ann2$V3))
mit <- data.frame(as.matrix(t2@assays$RNA@data[,pos]))

t2.ann2 <- t2.ann2[order(t2.ann2$V1),]
mit <- mit[,as.character(t2.ann2$V3)]
pos <- which(rownames(mat) %in% rownames(mit))
mat1 <- mat[pos,]
mit1 <- mit[rownames(mat1),]
mat2 <- cbind(mat1,mit1)
ann <- data.frame(clus = c(rep("T1",ncol(mat1)),as.character(t2.ann2$V1)))
rownames(ann) <- colnames(mat2)

pos <- which(t1.ann1$V2 %in% "CTCAGTCAGACAGTGTGACTGTGTCACA")
t1.ann4 <- t1.ann1[pos,]
pos <- which(t2.ann$V2 %in% "CTCAGTCAGACAGTGTGACTGTGTCACA") #CTGTGACACTCAGTGTGTCTGTCACACA GACTGTGACTCTGACAGTGAGTCTGTCT
t2.ann2 <- t2.ann[pos,]
table(t2.ann2$V1)

pos <- which(colnames(t1) %in% t1.ann4$V3)
mat <- data.frame(as.matrix(t1@assays$RNA@data[,pos]))
pos <- (which(colnames(t2) %in% t2.ann2$V3))
mit <- data.frame(as.matrix(t2@assays$RNA@data[,pos]))

t2.ann2 <- t2.ann2[order(t2.ann2$V1),]
mit <- mit[,as.character(t2.ann2$V3)]
pos <- which(rownames(mat) %in% rownames(mit))
mat1 <- mat[pos,]
mit1 <- mit[rownames(mat1),]
mat4 <- cbind(mat1,mit1)
ann2 <- data.frame(clus = c(rep("T1",ncol(mat1)),as.character(t2.ann2$V1)))
rownames(ann2) <- colnames(mat4)

pos <- which(rownames(mat4) %in% rownames(mat2))
mat5 <- mat4[pos,]
mat2 <- mat2[rownames(mat5),]
mat6 <- cbind(mat5,mat2)

ann2$clus <- paste0("barcode1_",ann2$clus)
ann$clus <- paste0("barcode2_",ann$clus)
ann3 <- rbind(ann2,ann)

clus <- unique(ann3$clus)
mitt <- matrix(nrow = nrow(mat6),ncol = length(clus))
for(i in 1:length(clus)){
  pos <- which(ann3$clus %in% clus[i])
  if(length(pos) > 1){
    matt <- mat6[,pos]
    mitt[,i] <- rowMeans(matt)
  }else{
    mitt[,i] <- as.numeric(as.character(mat6[,pos]))
  }
}
colnames(mitt) <- clus
pheatmap::pheatmap(cor(mitt),show_rownames = T,cluster_rows = T,cluster_cols = T)


################
#t1 <- FindClusters(t1,resolution = 0.2)
#t2 <- FindClusters(t2,resolution = 0.2)

pos <- which(t1$seurat_clusters %in% c(0,1,2))
cellnames <- colnames(t1)[pos]
t1 <- subset(t1,cells = cellnames)

t1.markers <- FindAllMarkers(t1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pos <- which(t1.markers$p_val_adj < 0.05)
t1.markers <- t1.markers[pos,]
table(t1.markers$cluster)
t1.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(t1, features = top10$gene) + NoLegend()

t2.markers <- FindAllMarkers(t2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pos <- which(t2.markers$p_val_adj < 0.05)
t2.markers <- t2.markers[pos,]
table(t2.markers$cluster)
t2.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(t2, features = top10$gene) + NoLegend()


clus <- unique(as.character(t1.markers$cluster))
feas <- unique(as.character(t2.markers$cluster))
mit <- matrix(nrow = length(clus),ncol = length(feas))
for(i in 1:length(clus)){
  for(j in 1:length(feas)){
    pos <- which(t1.markers$cluster %in% clus[i])
    t1.markers1 <- t1.markers[pos,]
    pos <- which(t2.markers$cluster %in% feas[j])
    t2.markers1 <- t2.markers[pos,]
    pos <- which(t1.markers1$gene %in% t2.markers1$gene)
    mit[i,j] <- length(pos) / nrow(t1.markers1) / nrow(t2.markers1)*1000
  }
}
rownames(mit) <- paste0("T1_",clus)
colnames(mit) <- paste0("T2_",feas)
pheatmap::pheatmap(mit,display_numbers = T)







