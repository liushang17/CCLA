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
pca <- readRDS("/Users/shangliu/01.terms/03.ecDNA/NC/01.Figure1/LSI.rds")
ob <- readRDS("/Users/shangliu/01.terms/03.ecDNA/NC/01.Figure1/ATAC.metadata.rds")

i <- "MDAMB453"
ITH <- c()
for (i in unique(ob$biostates)) {
temp_mat <- pca[which(ob$biostates==i),]
cutoff1 <- apply(temp_mat, 2, function(x) mean(x)-3*sd(x))
cutoff2 <- apply(temp_mat, 2, function(x) mean(x)+3*sd(x))
cell_bool <- apply(temp_mat, 1, function(x) sum(x[1:3] < cutoff1[1:3]) < 3 & sum(x[1:3] > cutoff2[1:3]) < 3)
temp_mat <- temp_mat[cell_bool,]
mean_m<-apply(temp_mat, 2, function(y) mean(y))
diversity <- mean(apply(temp_mat, 1, function(x) sqrt(sum((x-mean_m)^2))))
diversity
ITH <- c(ITH,diversity)

}
info <- data.frame(cellline = unique(ob$biostates),value = ITH)
info <- read.table(file = "/Users/shangliu/01.terms/03.ecDNA/NC/01.Figure1/ATAC.diversity.score.xls",sep = "\t",header = T)

info <- info[order(info$value),]
clus <- as.character(info$cellline)
library(ggplot2)
ggplot(info,aes(x=factor(cellline,levels = clus),y = value,fill = pattern))+geom_bar(stat ="identity")+theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+labs(x="",y="diversity score")
pos <- which(info$pattern %in% "Differential")
info1 <- info[pos,]
info2 <- info[-pos,]
t.test(info1$value,info2$value,alternative = "greater")

ggplot(info,aes(x = pattern,y = value,fill = pattern))+geom_violin()+geom_boxplot()+theme_classic()+labs(x="",y="diversity score")



