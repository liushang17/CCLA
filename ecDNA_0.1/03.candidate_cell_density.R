args<-commandArgs(T)
#library(ggpubr)
library(tidyr)
#cutoff=as.numeric(args[3])
#cutoff=0.05
cutoff=15
mydata<-read.table(args[1],header=F)
mydata_wide<-spread(mydata,V2,V3)
mydata_wide[is.na(mydata_wide)]<-0
rownames(mydata_wide)=mydata_wide[,1]
mydata_wide$V1<-NULL
mydata_wide[mydata_wide != 0 ]<-1
predensity<-data.frame(bin=rownames(mydata_wide))
predensity$cellnum<-apply(mydata_wide,1,sum)
#p1 <- ggdensity(predensity,x="cellnum",fill = "#c7e9c0",add = "mean", rug = TRUE)
#png(paste0(args[2],"_candidate_bin.png"),width=1200,height=600)
#print(p1)
#dev.off()
#ttt<-predensity[predensity$cellnum > cutoff*length(colnames(mydata_wide)),]
#ttt<-predensity[predensity$cellnum > cutoff*length(colnames(mydata_wide))& predensity$cellnum < (1-cutoff)*length(colnames(mydata_wide)),]
ttt<-predensity[predensity$cellnum > cutoff,]
ttt$ratio<-round(ttt$cellnum/length(colnames(mydata_wide)),2)
write.table(file= paste0(args[2],".candidate_bin.txt"),ttt,row.names=FALSE,col.names=FALSE,quote = FALSE,sep="\t")
