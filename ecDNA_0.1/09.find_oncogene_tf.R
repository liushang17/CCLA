args<-commandArgs(T)
.libPaths("/jdfssz1/ST_SUPERCELLS/P18Z10200N0350/USER/yuanyue/miniconda3/envs/r403/lib/R/library")
library(Seurat)
library(ggplot2)
library(RColorBrewer)

data<-read.table(args[1],stringsAsFactors=FALSE)
genelist<-data$gene_id
oncogene<-read.table("/hwfssz1/ST_SUPERCELLS/P21Z10200N0134/USER/zhaoxin/P18Z10200N0204/data_2/ATAC_g38/BGI500/04.ecDNA_detection/ref/Census_allTue.txt",header=T,sep="\t")
oo<-intersect(genelist,oncogene$Gene.Symbol)
write.table(oo,file=paste0(args[2],"_ecDNA_in_oncogene.txt"),quote=F,sep="\t",row.names=FALSE,col.names=FALSE)
q=length(oo)-1
#oncogene
m=723
#protein_coding_gene
n=50610
k=length(genelist)
#p=0.01,fish_test=-2
log(phyper(q,m,n,k,lower.tail = F),10)
length(oo)



tfs<-read.table("/hwfssz1/ST_SUPERCELLS/P21Z10200N0134/USER/zhaoxin/P18Z10200N0204/data_2/ATAC_g38/BGI500/04.ecDNA_detection/ref/trrust_rawdata.human.tsv",stringsAsFactors=FALSE)
ecDNA_tfs<-intersect(unique(tfs$V1),genelist)
write.table(ecDNA_tfs,file=paste0(args[2],"_ecDNA_in_tfs.txt"),quote=F,sep="\t",row.names=FALSE,col.names=FALSE)

ctl1<-readRDS(args[3])
ctl1.markers <- FindAllMarkers(ctl1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ctl1.markers<-ctl1.markers[which(ctl1.markers$avg_log2FC >0.5& ctl1.markers$p_val_adj<0.05& ctl1.markers$cluster=="Malignant"),]
pointgene<-intersect(ctl1.markers$gene,genelist)
write.table(pointgene,file=paste0(args[2],"_ecDNA_in_Malignant_markers.txt"),quote=F,sep="\t",row.names=FALSE,col.names=FALSE)

