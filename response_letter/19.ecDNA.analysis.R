ann <- read.table("/Users/shangliu/01.terms/03.ecDNA/NC/01.Figure1/table9.csv",sep = ",",header = T)
celllines <- unique(as.character(ann$Cell.Line))
mit <- matrix(nrow = length(celllines),ncol = 2)

for(i in 1:length(celllines)){
  pos <- which(ann$Cell.Line %in% celllines[i])
  ann1 <- ann[pos,]
  pos1 <- grep("NA",ann1$ecDNA.region)
  ann2 <- ann1[-pos1,]
  mit[i,1] <- celllines[i]
  mit[i,2] <- length(pos1) + length(unique(ann2$ecDNA.region))
}

celllins <- c("T-47D","SW620","SW480","SNB75","SF295","SF268","SCC4","RKO","MD-AMB-453","MD-AMB-361","MD-AMB-231","K-562",
              "Huh7","HS578T","HNSCCUM-02T","HNSCCUM-03T","HCT-8","HCT-15","HCT116","HCC1937","Hap1","Fadu","DLD1","COLO 205",
              "Caco2","BT-474")
celllins <- rev(celllins)
pos <- which(mit[,1] %in% celllins)
mit <- mit[pos,]

mit1 <- data.frame(mit)
rownames(mit1) <- mit1$X1
library(ggplot2)
ggplot(mit1,aes(x=factor(mit1$X1,levels = celllins),y = as.numeric(as.character(mit1$X2))))+geom_point()+theme_classic()+
  coord_flip()+labs(x="",y="")


library(ggplot2)
library(ggsignif)
setwd("/Users/shangliu/01.terms/03.ecDNA/NC/08.Siggenes")
file <- read.table("all_celllines_oncogene_varibility.txt")


ggplot(file,aes(V3,V4,fill=V3))+
  geom_violin(trim = F)+
  geom_boxplot(width=0.2,position=position_dodge(0.9),outlier.colour = NA,fill="white")+
  theme_bw()+
  theme(legend.position="none",
        axis.text.x = element_text(hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(hjust = 0.5, vjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text( size=rel(1)),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black",size=1))+
  scale_fill_manual(values = c("#d0d1e6","#0570b0"))+
  #scale_y_continuous(limits = c(15,30))+
  geom_signif(comparisons = list(c("oncogene","genes")),
              map_signif_level=T,
              textsize=6,test=t.test,step_increase=0.2)

all <- dir("./",pattern="*.oncogene_varibility.txt")
all1 <-all[-c(1,2,3,13,14)]
a <- list();for(i in 1:39){b <- read.table(paste0("/Users/shangliu/01.terms/03.ecDNA/NC/08.Siggenes/",all1[i]));a[[i]] <- b}
aa <- do.call(rbind,a)
#pdf("/jdfssz1/ST_SUPERCELLS/P21Z10200N0134/USER/zhangyuanhang1/CRC/CCLA/update_0425/supp/ecDNA_fig5/5C_update.pdf",width=4, height = 4)
ggplot(aa1,aes(V3,V4,fill=V3))+
  geom_violin(trim = F)+
  geom_boxplot(width=0.2,position=position_dodge(0.9),outlier.colour = NA,fill="white")+
  theme_bw()+
  theme(legend.position="none",
        axis.text.x = element_text(hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(hjust = 0.5, vjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text( size=rel(1)),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black",size=1))+
  scale_fill_manual(values = c("#d0d1e6","#0570b0"))+
  #scale_y_continuous(limits = c(15,30))+
  geom_signif(comparisons = list(c("oncogene","genes")),
              map_signif_level=F,
              textsize=6,test=t.test,step_increase=0.2)
#dev.off()

########
clus <- unique(as.character(ann$ecDNA.region))
pos <- grep("NA",clus)
clus <- clus[-pos]

ecdna <- NULL
for(i in 1:length(clus)){
  pos <- which(ann$ecDNA.region %in% clus[i])
  ann1 <- ann[pos,]
  pos <- which(aa$V2 %in% ann1$ecDNA.sub.region)
  aa1 <- aa[pos,]
  pos <- which(aa1$V4 == max(aa1$V4))
  cellt <- aa1$V2[-pos]
  ecdna <- c(ecdna,cellt)
}
pos <- which(aa$V2 %in% ecdna)
aa1 <- aa[-pos,]



preplot <- read.table("/Users/shangliu/01.terms/03.ecDNA/NC/lost/ecDNA_RNA_correlation_housegene_normalize_plotdata_exHKcellline.csv",sep = ",",header = T)
preplot$log2coverage <- log2(preplot$ecDNA_normalied + 1)
preplot$log2gene <- log2(preplot$RNA_normalied + 1)
library(ggpubr)
ggscatter(preplot,x = "log2coverage", y = "log2gene",add = "reg.line",conf.int = TRUE,
          size=1.5,add.params = list(color = "#854e4b",fill = "lightgray"))+
  stat_cor(method = "spearman")+
  theme_classic()+
  labs(x="ecDNA",y="Gene Expression")

pos <- grep("NA",ann$ecDNA.region)
ann1 <- ann[pos,]
ann2 <- ann[-pos,]
pos <- which(preplot$region %in% ann1$ecDNA.sub.region)
preplot1 <- preplot[pos,]
preplot2 <- preplot[-pos,]

mit <- NULL
clus <- unique(as.character(ann2$ecDNA.region))
for(i in 1:length(clus)){
  pos <- which(ann2$ecDNA.region %in% clus[i])
  ann3 <- ann2[pos,]
  pos <- which(preplot2$region %in% ann3$ecDNA.sub.region)
  if(length(pos) > 0){
    if(length(pos) == 1){
      tmp <- preplot2[pos,]
      mit <- rbind(mit,tmp)
    }else{
      tmp <- preplot2[pos,]
      tmp1 <- data.frame(region = clus[i], ecDNA_normalied = mean(tmp$ecDNA_normalied), RNA_normalied= mean(tmp$RNA_normalied),
                         sample = tmp$sample[i],log2coverage = 0,log2gene= 0)
      mit <- rbind(mit,tmp1)
  }
}
  
}
preplot3 <- rbind(preplot1,mit)

preplot3$log2coverage <- (log2(preplot3$ecDNA_normalied))
preplot3$log2gene <- (log2(preplot3$RNA_normalied))
library(ggpubr)
ggscatter(preplot3,x = "log2coverage", y = "log2gene",add = "reg.line",conf.int = TRUE,
          size=1.5,add.params = list(color = "#854e4b",fill = "lightgray"))+
  stat_cor(method = "spearman")+
  theme_classic()+
  labs(x="ecDNA",y="Gene Expression")



