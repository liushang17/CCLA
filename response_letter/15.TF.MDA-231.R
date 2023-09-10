library(Seurat)
library(ggsci)
library(ggrepel)
library(dplyr)
library(pheatmap)
library("reshape2")
library("tidyr")
library(data.table)
library(Matrix)
col <- c(pal_npg("nrc")(10)[1:6],"#FAFD7CFF","#FF6F00FF",pal_lancet("lanonc")(9)[c(1,3,7)],"#660099FF","#B5CF6BFF","#B24745FF","#CCFF00FF",
         "#FFCD00FF","#800000FF","#20854EFF","#616530FF","#FF410DFF","#EE4C97FF","#FF1463FF","#00FF00FF","#990080FF","#00FFFFFF",
         "#666666FF","#CC33FFFF","#00D68FFF","#4775FFFF","#C5B0D5FF","#FDAE6BFF","#79CC3DFF","#996600FF","#FFCCCCFF","#0000CCFF",
         "#7A65A5FF","#1A5354FF","#24325FFF")

load("/Users/shangliu/01.terms/03.ecDNA/NC/01.Figure1/allcellines_update_newbatch_filter_anno.Rdata")

ann <- ob@meta.data
sui <- data.frame(table(ann$Cell_line,ann$Lane))
pos <- which(sui$Freq == 0)
sui1 <- sui[-pos,]

ann$Cell_line <- gsub("_.*","",ann$Cell_line)
pos <- which(ann$Cell_line %in% "MDAMB231")
mat2 <- ob@assays$RNA@counts[,pos]

pbmc <- CreateSeuratObject(mat2)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features =   VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:30)

DimPlot(pbmc,label = T,cols = col)
VlnPlot(pbmc,features = "KLK3")

##################
genes <- read.table("/Users/shangliu/01.terms/03.ecDNA/NC/05.getgene/Nomal_MDAMB231_gene.txt",sep = "\t")
ase3 <- read.table("/Users/shangliu/01.terms/03.ecDNA/NC/08.ecDNA_cell_mtx/MDAMB231.ecDNA_cell_mtx.txt",sep = "\t",header = F)
tfinfo <- readRDS("/Users/shangliu/01.terms/03.ecDNA/NC/02.SNB75/3.1_regulons_forAUCell.Rds")
pos <- which(genes$V2 %in% "SOX10")
genes[pos,]
pos <- which(genes$V1 %in% "chr7:148300000_149900000")
genename <- genes[pos,2]
sed <- list(tes = genename)

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

meta <- read.table("/Users/shangliu/01.terms/03.ecDNA/NC/02.SNB75/MDAMB231.meta.xls",sep = "\t",header = T)
ggplot(meta,aes(x=UMAP_1,y=UMAP_2,color = Clusters2))+geom_point()
meta$cellname <- gsub(".*\\.","",rownames(meta))
pos <- which(meta$cellname %in% colnames(mat1))
meta1 <- meta[pos,]
mat2 <- mat1[,as.character(meta1$cellname)]

pos <- which(rownames(mat2) %in% "chr7:148300000_149900000")
meta1$value <- as.numeric(as.matrix(mat2[pos,]))
ggplot(meta1,aes(x=UMAP_1,y=UMAP_2,color = value))+geom_point()+ scale_color_gradient2(low = "lightblue",mid = "white",high = "red")

ggplot(meta1,aes(x=as.character(Clusters2),y = value))+geom_violin()+geom_boxplot()+theme_classic()

pos <- which(meta1$Clusters2 %in% "0")
meta2 <- meta1[pos,]
meta3 <- meta1[-pos,]
t.test(meta2$value,meta3$value)

meta1$predictedCell_Un <- (paste0("MDAMB231_2_",meta1$predictedCell_Un))
cellnames <- unique(unique(meta1$predictedCell_Un))

pbmc1 <- subset(pbmc,cells = cellnames)
DimPlot(pbmc1,label = T,cols = col)
meta2 <- unique(meta1[,c("predictedCell_Un","Clusters2")])
rownames(meta2) <- 1:nrow(meta2)
pbmc1$new_type <- "R0"
for(i in 1:length(colnames(pbmc1))){
  pos <- which(meta2$predictedCell_Un %in% colnames(pbmc1)[i])
  pbmc1$new_type[i] <- meta2$Clusters2[pos[1]]
}
DimPlot(pbmc1,label = T,cols = col,group.by = "new_type")
pbmc1 <- AddModuleScore(pbmc1,features = sed)
VlnPlot(pbmc1,features = "Cluster1",group.by = "new_type",cols = rep("white",5),pt.size = 0)+geom_boxplot()
#clus.mak <- FindMarkers(pbmc1,group.by = "new_type",ident.1 = "0")

pos <- which(pbmc1$new_type %in% "0")
tmp <- pbmc1$Cluster1[pos]
pos <- which(pbmc1$new_type %in% "1")
tmp1 <- pbmc1$Cluster1[pos]
t.test(tmp,tmp1,alternative = "greater")

pos <- which(rownames(pbmc1@assays$RNA@data) %in% "TP53")
tmpgene <- as.numeric(as.character(as.matrix(pbmc1@assays$RNA@data[pos,])))
pos <- which(pbmc1$new_type %in% "0")
tmp <- tmpgene[pos]
pos <- which(pbmc1$new_type %in% "2")
tmp1 <- tmpgene[pos]
t.test(tmp,tmp1,alternative = "greater")

pos <- which(tfinfo$TF %in% "SOX10")
tfinfo[pos,]

tes <- readRDS("/Users/shangliu/01.terms/03.ecDNA/NC/02.SNB75/3.4_regulonAUC.MD231.Rds")
tes1 <- data.frame(t(as.matrix(tes@assays@data$AUC)))
rownames(tes1) <- paste0("MDAMB231_2_",rownames(tes1))
tes2 <- tes1[colnames(pbmc1),]
pbmc1@meta.data <- cbind(pbmc1@meta.data,tes2)
sed <- list(tes = as.character(strsplit("CASP3/FOXO3/EZH2/RECQL4/FAS/TCF7L2/CDKN1B/CDK4/AKT1/FUS/TP53/STAT3/BRCA1/BCL3/BAX/POLD1/CHEK2/ELF4",split = "/")[[1]]))
