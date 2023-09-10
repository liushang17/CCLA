library(ggplot2)
library(org.Hs.eg.db)
library(clusterProfiler)
hallmarks <- read.gmt("/hwfssz1/ST_SUPERCELLS/P21Z10200N0134/USER/zhaoxin3/P18Z10200N0204/data_2/cellline/marker/result/result2/h.all.v7.2.symbols.gmt")

# MDAMB231 FOXA2-----------------
regulons <-readRDS("/mydata/MDAMB231/SCENIC/int/2.6_regulons_asGeneSet.Rds")
M_c0<-regulons[c("FOXA2")]
M_c0
ego_ALL<- enricher(M_c0$FOXA2, TERM2GENE = hallmarks)
ego_ALL@result<-ego_ALL@result[ego_ALL@result$pvalue<=0.05,]
pdf("FOXA2_dotplot_pvalue.pdf",5,3)
ggplot(ego_ALL@result, aes(x=GeneRatio, y=reorder(ID,Count))) +geom_point(aes(color=pvalue, size=Count))+ylab("")
dev.off()

# RPMI8226 NFE2L2
data <- read.table("/jdfssz1/ST_SUPERCELLS/P21Z10200N0134/USER/zhangyuanhang1/CRC/CCLA/update_0425/01.result/umap.123/20230322/test/pyscenic/RPMI8226/RPMI8226_regulons.tsv",
	header=T,sep="\t")
data1 <- data[grepl("NFE2L2",data$X),]
gene <- as.list(data1[,9])
gene1 <- sapply(gene,function(x){
  tmp1 <- gsub(",| |0.[0-9]+","",x)
  tmp2 <- strsplit(tmp1,split="\\)\\(")
  tmp2 <- gsub("\\[\\(|\\)\\]","",tmp2[[1]])


})
gene2 <- Reduce(union,gene1)
ego_ALL<- enricher(gene2, TERM2GENE = hallmarks)
ego_ALL@result<-ego_ALL@result[ego_ALL@result$pvalue<=0.05,]

pdf("NFE2L2_dotplot_pvalue.pdf",5,3)
ggplot(ego_ALL@result, aes(x=GeneRatio, y=reorder(ID,Count))) +geom_point(aes(color=pvalue, size=Count))+ylab("")
dev.off()
# SNB75 NFYB ------------------
data <- read.table("/jdfssz1/ST_SUPERCELLS/P21Z10200N0134/USER/zhangyuanhang1/CRC/CCLA/update_0425/01.result/umap.123/20230322/test/pyscenic/SNB75/SNB75_regulons.tsv",
	header=T,sep="\t")
data1 <- data[grepl("NFYB",data$X),]
gene <- as.list(data1[,9])
gene1 <- sapply(gene,function(x){
  tmp1 <- gsub(",| |0.[0-9]+","",x)
  tmp2 <- strsplit(tmp1,split="\\)\\(")
  tmp2 <- gsub("\\[\\(|\\)\\]","",tmp2[[1]])


})
gene2 <- Reduce(union,gene1)
ego_ALL<- enricher(gene2, TERM2GENE = hallmarks)
ego_ALL@result<-ego_ALL@result[ego_ALL@result$pvalue<=0.05,]
pdf("NFYB_dotplot_pvalue.pdf",8,4)
ggplot(ego_ALL@result, aes(x=GeneRatio, y=reorder(ID,Count))) +geom_point(aes(color=pvalue, size=Count))+ylab("")
dev.off()