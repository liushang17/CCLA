args<-commandArgs(T)
#library(ggplot2)
#library(tidyr)
data<-read.table(args[1],header=F)
candidate=data
#cell_qc<-read.table(args[2],header=T,row.names=1,sep=",")
gene<-read.table(args[2],header=F)
#cell_qc<-cell_qc[which(rownames(cell_qc)%in%candidate$V2),]
#cell_qc_tss<-data.frame(region="whole",cell=rownames(cell_qc),tssproportion=cell_qc$tssProportion)
tss<-read.table("/hwfssz1/ST_SUPERCELLS/P21Z10200N0134/USER/zhaoxin/P18Z10200N0204/data_2/ATAC_g38/BGI500/02.cellline/No_Hypoxia/test_ecDNA/bin_test/bin/tss.bed",header=F)
#tss <- read.table("/hwfssz5/ST_SUPERCELLS/P21Z10200N0134/19.ecDNA/02.method/ref/mm10/tss.bed",header=F)
tss<-tss[tss$V7== "protein_coding",]
names(candidate)<-c("region","cell","tss","sum","tssproportion")
#cellcut<-length(unique(candidate$cell))*0.05
cellcut<-15
candidate<-candidate[which(candidate$sum > 8),]
candidate<-candidate[which(candidate$tssproportion <= 0.25),]

tt<-data.frame(rowSums(table(candidate$region,candidate$cell)))
names(tt)<-c("cellnumber")
tt$region<-rownames(tt)
tt<-tt[which(tt$cellnumber>cellcut),]

candidate_2<-candidate[which(candidate$region %in% tt$region),]

data<-data[which(data$V1 %in% tt$region),]
ecincell<-data.frame(region=candidate_2$region,cell=candidate_2$cell,judge="1")
mylist<-c(paste0(candidate_2$region,"_",candidate_2$cell))

for (i in 1:nrow(data)){
	t<-paste0(data[i,]$V1,"_",data[i,]$V2)
	if (!t %in% mylist) {
	temp<-data.frame(region=data[i,]$V1,cell=data[i,]$V2,judge="0")
	ecincell<-rbind(ecincell,temp)
	}
}
write.table(ecincell,paste0(args[3],".region_cell_judge.txt"),quote=F,sep="\t",row.names=F)


tt$region<-NULL
#candidate$tss<-NULL
#candidate$sum<-NULL
#mydata<-rbind(candidate,cell_qc_tss)
#mydata<-candidate
#test_region<-sample(levels(mydata$region), 15, replace = FALSE)
#test_region[16]<-"whole"
#sample1<-mydata[which(mydata$region %in% test_region),]

#cellm<-tapply(mydata$tssproportion,mydata$region,mean)
#cellmd<-tapply(mydata$tssproportion,mydata$region,median)
#cellm<-as.data.frame(cellm)
#cellmd<-as.data.frame(cellmd)
#savedata<-cbind(cellm,cellmd)
#savedata<-savedata[savedata$cellm<0.1,]
write.table(tt,paste0(args[3],".filterd_region.txt"),quote=F,sep="\t")

names(gene)<-c("region","gene_id")
gene<-gene[gene$region %in% rownames(tt),]
gene<-gene[gene$gene_id %in% tss$V4,]
write.table(gene,paste0(args[3],".final_codinggene.txt"),quote=F,sep="\t")
#names(cellm)<-c("mean")
#cellm$region<-rownames(cellm)


#pdf(paste0(args[3],"_box_line.pdf"),width=,height=)
#ggplot(sample1, aes(x=region,y=tssproportion)) +geom_boxplot(aes(fill=region))+labs(x="Regions", y = "TSS proportion", fill = "Regions")+theme(plot.title = element_text(hjust = 0.8),axis.text.x=element_text(angle=45),legend.position = "bottom")
#ggplot(cellm,aes(x=region,y=mean,group=1))+geom_line()
#dev.off()
