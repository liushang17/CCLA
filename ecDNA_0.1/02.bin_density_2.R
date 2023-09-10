args<-commandArgs(T)
#library(ggpubr)
library(tidyr)
library(pheatmap)
mydata<-read.table(args[1],header=F)
title <-args[2]
#p1 <- ggdensity(mydata,x="V3",fill = "lightgray",add = "mean", rug = TRUE)
#pdf(args[2],height=6,width=10)
#print(p1)
#dev.off()
#setwd("/hwfssz1/ST_PRECISION/USER/zhaoxin/P18Z10200N0204/data_2/ATAC_g38/BGI500/02.cellline/No_Hypoxia/test_ecDNA/bin_test/bin/test/")
#mydata$V4<-exp(mydata$V3)
test_scCell<-sample(mydata$V2, 10, replace = FALSE)
my_newdata<-mydata[which(mydata$V2 %in% test_scCell),]
#p2 <- ggdensity(my_newdata,x="V3",color="V2", rug = TRUE)
#png(paste0(title,"_density_test",".png"),width=1200,height=600)
#print(p2)
#dev.off()

mydata$V4<-mydata$V3
mydata$V3<-NULL
result<-NULL
preheatmap_result<-NULL
for (cell in levels(factor(mydata$V2))){
	my_newdata<-mydata[mydata$V2 == cell,]
	MaxY1_index <- which.max(density(my_newdata$V4)$y)
	MaxX1 <- density(my_newdata$V4)$x[MaxY1_index]
	candidate_df<-my_newdata[my_newdata$V4/MaxX1 >= 5.5,]
	#my_newdata[which(my_newdata$V4/MaxX1 < 2.5),]$V4 = 0
	result<-rbind(result,candidate_df)
	#preheatmap_result<-rbind(preheatmap_result,my_newdata)
}

write.table(file = paste0(title,".candidate.txt"),result,row.names=FALSE,col.names=FALSE,quote = FALSE,sep="\t")

