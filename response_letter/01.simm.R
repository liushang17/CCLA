library(Seurat)
library(ggsci)
library(ggrepel)
library(dplyr)
col <- c(pal_npg("nrc")(10)[1:6],"#FAFD7CFF","#FF6F00FF",pal_lancet("lanonc")(9)[c(1,3,7)],"#660099FF","#B5CF6BFF","#B24745FF","#CCFF00FF",
         "#FFCD00FF","#800000FF","#20854EFF","#616530FF","#FF410DFF","#EE4C97FF","#FF1463FF","#00FF00FF","#990080FF","#00FFFFFF",
         "#666666FF","#CC33FFFF","#00D68FFF","#4775FFFF","#C5B0D5FF","#FDAE6BFF","#79CC3DFF","#996600FF","#FFCCCCFF","#0000CCFF",
         "#7A65A5FF","#1A5354FF","#24325FFF")
load("/Users/shangliu/01.terms/03.ecDNA/NC/01.Figure1/allcellines_update_newbatch_filter_anno.Rdata")

ann <- ob@meta.data
sui <- data.frame(table(ann$Cell_line,ann$Lane))
pos <- which(sui$Freq == 0)
sui1 <- sui[-pos,]

metadata <- read.table("/Users/shangliu/01.terms/03.ecDNA/NC/01.Figure1/Metadata.xls",sep = "\t",header = T)
bulk <- readRDS("/Users/shangliu/01.terms/03.ecDNA/NC/01.Figure1/OmicsExpressionGenesExpectedCountProfile.rds")

pos <- which(sui1$Var1 %in% metadata$StrippedCellLineName)
sui2 <- sui1[pos,]
sui3 <- sui1[-pos,]

run <- "CL200136256_L01"
pos <- which(sui2$Var2 %in% run) 
sui3 <- sui2[pos,]

pos <- which(metadata$StrippedCellLineName %in% sui3$Var1)
metadata1 <- metadata[pos,]
pos <- which(rownames(bulk) %in% metadata1$ProfileID)
bulk1 <- data.frame(t(bulk[pos,]))
bulk1$genename <- gsub("..EN.*","",rownames(bulk1))
ces <- data.frame(table(bulk1$genename))
pos <- which(ces$Freq > 1)
ces1 <- ces[-pos,]
pos <- which(bulk1$genename %in% ces1$Var1)
bulk2 <- bulk1[pos,]
rownames(bulk2) <- bulk2$genename


pos <- which(ob$Cell_line %in% sui3$Var1)
cellname <- colnames(ob)[pos]
ob1 <- subset(ob,cells = cellname)
mat <- as.matrix(ob1@assays$RNA@counts)
pos <- which(rownames(mat) %in% rownames(bulk2))
mat1 <- mat[pos,]
bulk2 <- bulk2[rownames(mat1),]
bulk2 <- bulk2[,-c(4)]
mat2 <- cbind(mat1,bulk2)

pbmc <- CreateSeuratObject(mat2)
pbmc$cellname <- c(as.character(ob1$Cell_line),"HCT116","MCF7","A253")
#pbmc$cellname <- c(as.character(ob1$Cell_line),"K562","SW620","MDAMB361")
pbmc$type <- c(rep("single",ncol(mat1)),"bulk","bulk","bulk")

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features =   VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)

DimPlot(pbmc,label = T,cols = col)
DimPlot(pbmc,label = T,group.by = "cellname",cols = col)
DimPlot(pbmc,label = T,group.by = "type",cols = col)
############
bulk_exp <- as.matrix(pbmc@assays$RNA@data[,1227:1229])
pbmc1 <- subset(pbmc,cells = colnames(pbmc)[1:1226])
#bulk_exp <- as.matrix(pbmc@assays$RNA@data[,1335:1337])
#pbmc1 <- subset(pbmc,cells = colnames(pbmc)[1:1334])

Idents(pbmc1) <- pbmc1$cellname
ob1.ave <- AverageExpression(pbmc1,return.seurat = T)
exp2 <- ob1.ave@assays$RNA@data
exp2 <- exp2[rownames(bulk_exp),]
exp3 <- cbind(exp2,bulk_exp)
pheatmap::pheatmap(exp3,cluster_rows = F,show_rownames = F,scale = "column")
genename <- pbmc@assays$RNA@var.features
pheatmap::pheatmap(exp3[genename,],cluster_rows = F,show_rownames = F,scale = "column")
suicor <- cor(exp3[genename,])
suicor <- suicor[1:3,4:6]
pheatmap::pheatmap(suicor)

Idents(pbmc) <- pbmc$cellname
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
exp4 <- exp3[top10$gene,]
na <- colnames(exp4)
na <- na[c(1,6,4,2,3,5)]
pheatmap::pheatmap(exp4[,na],cluster_rows = F,cluster_cols = F)

library(clusterProfiler)
genes <- read.gmt("/Users/shangliu/01.terms/03.ecDNA/NC/01.Figure1/gene_set_library_up_crisp.gmt")
pos <- which(genes$term %in% c("HCT116","MCF7","A253"))
genes1 <- genes[pos,]
pos <- which(genes1$gene %in% rownames(exp3))
genes1 <- genes1[pos,]
exp4 <- exp3[as.character(genes1$gene),]
na <- colnames(exp4)
na <- na[c(1,6,4,2,3,5)]
pheatmap::pheatmap(exp4[,na],cluster_rows = F,cluster_cols = F)

pos <- which(pbmc.markers$gene %in% genes1$gene)
pbmc.markers1 <- pbmc.markers[pos,]
pbmc.markers1 %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
exp4 <- exp3[top10$gene,]
na <- colnames(exp4)
na <- na[c(1,6,4,2,3,5)]
pheatmap::pheatmap(exp4[,na],cluster_rows = F,cluster_cols = F)

genename <- c("KHDC1L","CHURC1","PCDHB13","LINC00294","IMMP1L","JKAMP","FBXO3","ALOX12","DCAF4","ERH","PRRG4","ELP4",
              "OR51B2","ABHD10","SSFA2","DNAAF3","OR51I1","OR51I1","RAC3","VTI1A","MCMBP","PDP1",
              "SGCG","DSCAM","PARD6B","ELP2","RGS22","LINC00052","COX6C","DNASE1L2","SFXN1","GPATCH2")
pos <- which(genename %in% rownames(exp3))
genename <- genename[pos]
exp4 <- exp3[genename,]
na <- colnames(exp4)
na <- na[c(1,6,4,2,3,5)]
pheatmap::pheatmap(exp4[,na],cluster_rows = F,cluster_cols = F)

#############
pos <- grep("Breast|Neck|Colorectal",metadata$OncotreePrimaryDisease)
metadata2 <- metadata[-pos,]
pos <- which(rownames(bulk) %in% metadata2$ProfileID)
bulk1 <- data.frame(t(bulk[pos,]))
bulk1$genename <- gsub("..EN.*","",rownames(bulk1))
ces <- data.frame(table(bulk1$genename))
pos <- which(ces$Freq > 1)
ces1 <- ces[-pos,]
pos <- which(bulk1$genename %in% ces1$Var1)
bulk2 <- bulk1[pos,]
rownames(bulk2) <- bulk2$genename
bulk2 <- bulk2[rownames(mat1),]
bulk2 <- bulk2[,-c(1264)]
sed <- sample(1:ncol(bulk2),3)
bulk3 <- bulk2[,sed]
mat3 <- cbind(mat2,bulk3)

pbmc <- CreateSeuratObject(mat3)
pbmc$cellname <- c(as.character(ob1$Cell_line),"HCT116","MCF7","A253",rep("other",ncol(bulk3)))
pbmc$type <- c(rep("single",ncol(mat1)),"bulk","bulk","bulk",rep("other",ncol(bulk3)))

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features =   VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)


DimPlot(pbmc,label = T,cols = col)
DimPlot(pbmc,label = T,group.by = "cellname",cols = col)
DimPlot(pbmc,label = T,group.by = "type",cols = col)






