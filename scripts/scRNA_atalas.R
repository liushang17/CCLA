library(Seurat)
library(ggplot2)
load("mydata/allcellines_scRNA.Rdata")
# umap plot by cell line
anno_new <- read.csv("allcellines_color_cancer_aliasname.csv",row.names=1)# color 
Idents(ob) <- "Alias_update"
colorlable<-anno_new[levels(ob@active.ident),]$color
pdf("all_uamp_RNA_new.pdf",width = 10,height = 6)
DimPlot(ob, reduction = "umap",label = T,
        cols = colorlable)+
  #theme_classic(base_size = 4)+
  theme(axis.line = element_line(size = 1, colour = "black"))+
  guides(colour=guide_legend(override.aes=list(size=2)))+
  NoLegend()
dev.off()
#boxplot by UMI
data2 <- data.frame(Cell_lines=ob$Alias_update,nUMI=ob$nCount_RNA)
data2$Cell_lines <- factor(data2$Cell_lines,levels = levels(ob@active.ident))

data2$UMI <- log10(data2$nUMI)
pdf("nUMI_boxplot_update.pdf",width=8,height=5)
ggplot(data2, aes(x=Cell_lines, y=UMI, fill=Cell_lines)) + 
	geom_boxplot(outlier.colour = NA,color="#a39e9e")+
	theme_classic()+scale_fill_manual(values=colorlable)+
	ylab("Number of UMI (log10)")+
	theme(legend.position="none")+theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1))+scale_y_continuous(limits=c(2.5, 5.5))
dev.off()

# boxplot by gene 
data2 <- data.frame(Cell_lines=ob$Alias_update,nGene=ob$nFeature_RNA)
data2$Cell_lines <- factor(data2$Cell_lines,levels = levels(ob@active.ident))
data2$nGene <- log10(data2$nGene)
pdf("nGene_boxplot_update.pdf",width=8,height=5)
ggplot(data2, aes(x=Cell_lines, y=nGene, fill=Cell_lines))+  
	geom_boxplot(outlier.colour = NA,color="#a39e9e")+
	theme_classic()+scale_fill_manual(values=colorlable)+
	theme(legend.position="none")+ylab("Number of gene (log10)")+
	theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1))+scale_y_continuous(limits=c(2.5, 4.3))
dev.off()

# boxplot by cell number 
data <- as.data.frame(table(ob$Alias_update))
colnames(data) <- c("Cell_lines","nCell")
data$Cell_lines <- factor(data$Cell_lines,levels = levels(ob@active.ident))
pdf("nCell_boxplot_update.pdf",width=8,height=5)
ggplot(data, aes(x=Cell_lines, y=nCell, fill=Cell_lines)) + geom_bar(stat='identity') + 
	theme_classic()+scale_fill_manual(values=colorlable)+theme(legend.position="none")+
	ylab("Number of cell")+
	theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1))+scale_y_continuous(limits=c(0, 1500))
dev.off()

# umap plot by batch
Idents(ob)<-"cols_tmp"

cols <- c(Caco2 = "#5c4f74",Caco2_2 = "#4f745c",
        MDAMB231 = "#f9615f",MDAMB231_2 = "#615ff9",SCC4 = "#ffad5e",SCC4_2 = "#5effad",
        Others = "#d3d3d3")
pdf("batch_update.pdf",width = 8,height = 6)

DimPlot(ob, reduction = "umap",label = F,
        cols = cols)+
  theme(axis.line = element_line(size = 1, colour = "black"))+
  guides(colour=guide_legend(override.aes=list(size=2)))
dev.off()

# umap plot by nCount_RNA
pdf("UMI_umap_update.pdf",5,4)
ob$"log10(nCount_RNA)" <- log10(ob$nCount_RNA)
FeaturePlot(ob,features="log10.nCount_RNA.",cols=c("#ecd696","#914d43"),pt.size=0.5)+
dev.off()


# umap plot by library 
color <- sample(coloranno,21)
cols <- unique(ob$Library)
names(color) <- unique(ob$Library)
Idents(ob) <- "Library"
pdf("library_umap_update.pdf",8,6)
DimPlot(ob, reduction = "umap",label = F,
        cols = color)+
  theme(axis.line = element_line(size = 1, colour = "black"))+
  guides(colour=guide_legend(override.aes=list(size=2)))
 dev.off()
# umap plot  by cancer type 
Idents(ob) <- "Cancer_type"
color <- c("#67b988","#df6f73","#9ad2a7","#a983a8","#f2ea93","#d8a38e","#3c8293","#e53b3e","#f48565","#b5b6b9")  
names(color) <-unique(ob$Cancer_type)[order(unique(ob$Cancer_type))]

pdf("Cancer_type_umap_update.pdf",8,6)
DimPlot(ob, reduction = "umap",label = F,
        cols = color)+
  theme(axis.line = element_line(size = 1, colour = "black"))+
  guides(colour=guide_legend(override.aes=list(size=2)))
dev.off()

