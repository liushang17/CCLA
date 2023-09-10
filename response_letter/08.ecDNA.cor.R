library(pheatmap)
library("reshape2")
library("tidyr")
library(data.table)
library(Matrix)
library(ggplot2)
indir <- "/Users/shangliu/01.terms/03.ecDNA/NC/tmp/"
files <- list.files(indir)

#numdi <- matrix(nrow = length(files),ncol = 3)
#for(m in 1:length(files)){
files <- "Normal_CL82_HT29.ecDNA_cell_mtx.txt"
  filename <- paste0(indir,files)
  ase3 <- read.table(filename,sep = "\t")
  
  ###### trans
  
  ase3$FC <- ase3$V3
  ase3$cellID<-ase3$V2
  ase3$geneID <- ase3$V1
  pos <- which(ase3$symbol %in% NA)
  if(length(pos) > 0){ase3 <- ase3[-pos,]}
  gene=unique(ase3$geneID)
  cell=unique(ase3$cellID)
  gene_idx=c(1:length(gene))
  cell_idx=c(1:length(cell))
  names(gene_idx)=gene
  names(cell_idx)=cell
  ase3$num <- 1
  mat=sparseMatrix(i=gene_idx[ase3$geneID],j=cell_idx[ase3$cellID],x=ase3$FC)
  mit=sparseMatrix(i=gene_idx[ase3$geneID],j=cell_idx[ase3$cellID],x=ase3$num)
  mit[mit == 0] <- 1
  
  mat1 <- mat / mit
  mat1 <- as.matrix(mat1)
  mat1[mat1 == NA] <- 0
  
  rownames(mat1)=gene
  colnames(mat1)=cell
  if(nrow(mat1) > 5){
  #############
  mat2 <- cor(t(mat1),method = "pearson")
  p1 <- pheatmap::pheatmap(mat2,show_rownames = F,show_colnames = F)
  
  num <- ceiling(ncol(mat2)/5)
  clusmat <- data.frame(cutree(p1$tree_row,k = num))
  
  info <- matrix(nrow = num,ncol = 3)
  for(i in 1:num){
    pos <- which(clusmat$cutree.p1.tree_row..k...num. == i)
    tmpname <- rownames(clusmat)[pos]
    if(length(pos) > 1){
      pos1 <- which(rownames(mat2) %in% tmpname)
      matt <- mat2[pos1,pos1]
      pos2 <- which(matt > 0.3)
      info[i,1] <- i
      info[i,2] <- (length(pos2) - nrow(matt) ) /2
      info[i,3] <- (length(pos2) - nrow(matt)) / (nrow(matt) * (nrow(matt) - 1)) 
    }else{
      info[i,1] <- i
      info[i,2] <- 1
      info[i,3] <- 0
    }
  }
  
  #plot(density(info[,3]))+abline(v = 0.7)
 
  clusmat$new_clus <- clusmat$cutree.p1.tree_row..k...num.
  pos <- which(info[,3] > 0.7)
  pos <- which(clusmat$new_clus %in% info[pos,1])
  clusmat$new_clus[-pos] <- "Other"
  
  clusmat <- clusmat[order(clusmat$new_clus),]
  clusmat1 <- data.frame(clus = clusmat$new_clus)
  rownames(clusmat1) <- rownames(clusmat)
  
  mat3 <- mat2[rownames(clusmat1),rownames(clusmat1)]
  #pheatmap::pheatmap(mat3,annotation_row = clusmat1,cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F)
  
  #length(unique(clusmat1$clus))
  pos1 <- which(clusmat1$clus %in% "Other")
  numdi[m,1] <- gsub(".ecDNA_cell_mtx.txt","",files[m])
  numdi[m,2] <- length(unique(ase3$V1))
  numdi[m,3] <- length(pos1) + length(unique(clusmat1$clus)) - 1
  }else{
    numdi[m,1] <- gsub(".ecDNA_cell_mtx.txt","",files[m])
    numdi[m,2] <- nrow(mat1)
    numdi[m,3] <- nrow(mat1)
  }
  
#}

bo <- c(2,3,)

ase3 <- read.table("/Users/shangliu/01.terms/03.ecDNA/NC/01.Figure1/Nomal_MDAMB231.ecDNA_cell_mtx.txt",sep = "\t")


mat2 <- cor(t(mat1),method = "pearson")
p1 <- pheatmap::pheatmap(mat2,show_rownames = F,show_colnames = F)

num <- ceiling(ncol(mat2)/5)
clusmat <- data.frame(cutree(p1$tree_row,k = num))

info <- matrix(nrow = num,ncol = 3)
for(i in 1:num){
  pos <- which(clusmat$cutree.p1.tree_row..k...num. == i)
  tmpname <- rownames(clusmat)[pos]
  if(length(pos) > 1){
  pos1 <- which(rownames(mat2) %in% tmpname)
  matt <- mat2[pos1,pos1]
  pos2 <- which(matt > 0.3)
  info[i,1] <- i
  info[i,2] <- (length(pos2) - nrow(matt) ) /2
  info[i,3] <- (length(pos2) - nrow(matt)) / (nrow(matt) * (nrow(matt) - 1)) 
  }else{
    info[i,1] <- i
    info[i,2] <- 1
    info[i,3] <- 0
  }
}

plot(density(info[,3]))+abline(v = 0.7)

pos <- which(clusmat$cutree.p1.tree_row..k...num. == 88)
tmpname <- rownames(clusmat)[pos]
pos1 <- which(rownames(mat2) %in% tmpname)
matt <- mat2[pos1,pos1]
pheatmap::pheatmap(matt,show_rownames = T,show_colnames = F,display_numbers = T) 

clusmat$new_clus <- clusmat$cutree.p1.tree_row..k...num.
pos <- which(info[,3] > 0.7)
info1 <- info[pos,]
pos <- which(clusmat$new_clus %in% info1[,1])
clusmat$new_clus[-pos] <- "Other"

clusmat <- clusmat[order(clusmat$new_clus),]
clusmat1 <- data.frame(clus = clusmat$new_clus)
rownames(clusmat1) <- rownames(clusmat)

mat3 <- mat2[rownames(clusmat1),rownames(clusmat1)]
pheatmap::pheatmap(mat3,annotation_row = clusmat1,cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F)

length(unique(clusmat1$clus))
pos1 <- which(clusmat1$clus %in% "Other")
length(pos1) + length(unique(clusmat1$clus)) - 1

pos <- which(clusmat$cutree.p1.tree_row..k...num. == 88)
tmpname <- rownames(clusmat)[pos]
pos1 <- which(rownames(mat2) %in% tmpname)
matt <- mat2[pos1,pos1]
#pheatmap::pheatmap(matt,show_rownames = T,show_colnames = F,display_numbers = T) 

sui <- read.table("/Users/shangliu/01.terms/03.ecDNA/NC/02.SNB75/All.num.txt",sep = "\t",header = T)
sui1 <- data.frame(t(sui[,1:31]))
pos <- which(sui1$X2 > 0)
sui2 <- sui1[pos,]
sui3 <- data.frame(clus = rep(c("Before","After"),each = nrow(sui1)), value = c(sui1$X1,sui1$X2))
pos <- which(sui3$value > 0)
sui3 <- sui3[pos,]
ggplot(sui3,aes(x=factor(clus,levels = c("Before","After")),y = value))+geom_boxplot()+theme_classic()+labs(x="",y="")


