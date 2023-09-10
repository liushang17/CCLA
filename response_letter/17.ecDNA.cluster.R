mat <- read.table("/Users/shangliu/01.terms/03.ecDNA/NC/ecDNA.set.all2.xls",sep = "\t")
sui <- data.frame(table(mat$type))
pos <- grep("Hap1",sui$Var1)
sui[pos,]

cellline <- c("Normal_CL200149546_L01_HCT15","Normal_CL81_HCT8","Normal_CL82_HT29","Normal_CL200149516_L01_SW620","Normal_CL81_SW480",
              "Normal_Lovo","Normal_RKO","Nomal_MDAMB231","Normal_CL82_MDAMB361","Normal_CL200149546_L01_MDAMB453","Normal_CL82_MDAMB468",
              "Normal_CL16_BT474","Normal_CL200149516_L01_T47D","Normal_CL81_HCC1937","Fadu_CL100169139_L02_C7","SCC4_CL100169089_L01_C8",
              "K562_862","RPMI8226_CL100169086_L01_C4","SF268_CL100169137_L02_C6","SNB75_CL100169085_L02_C6","SF295_372","Huh7_2","786-0_CL100169139_L02_C3",
              "Hap1_CL100169139_L01_C9","Nomal_CL801_Colo205","Hunkel_CL100169089_L02_C3","Deucer_CL100169137_L01_C2","Normal_CL200149546_L01_DLD1","CaCo2","HCT116","HS578T")
cellname <- c("HCT-15","HCT-8","HT-29","SW620","SW480","LoVo","RKO","MD-AMB-231","MD-AMB-361","MD-AMB-453","MD-AMB-468","BT-474",
              "T-47D","HCC1937","Fadu","SCC4","K-562","RPMI 8226","SF268","SNB75","SF295","Huh7","786-O","Hap1","COLO 205","HNSCCUM-02T","HNSCCUM-03T","DLD1","Caco2","HCT116","HS578T")

pos <- which(sui$Var1 %in% cellline)
sui1 <- sui[-pos,]

pos <- which(mat$type %in% cellline)
mat1 <- mat[pos,]
mat1$cellname <- 0

for(i in 1:length(cellline)){
  pos <- which(mat1$type %in% cellline[i])
  mat1$cellname[pos] <- cellname[i]
}

mat2 <- mat1[,c(4,3,1)]
mat2$cluster <- paste0(mat2$cellname,":",mat2$clus)
#write.table(mat2,file = "/Users/shangliu/01.terms/03.ecDNA/NC/ecDNA.set.all.new.xls",sep = "\t",quote=F,row.names = F)
##########
mat.new <- NULL
mat2$chr <- gsub(":.*","",mat2$ecDNA)
for( i in 1:length(cellname)){
  pos <- which(mat2$cellname %in% cellname[i])
  mat3 <- mat2[pos,]
  clus <- unique(mat3$clus)
  pos1 <- which(clus %in% "other")
  pos2 <- which(clus %in% "Other")
  if(length(c(pos1,pos2)) > 0){
    clus <- clus[-c(pos1,pos2)]
    pos <- which(mat3$clus %in% c("other","Other"))
    mat4 <- mat3[pos,]
    mat.new <- rbind(mat.new,mat4)
  }
  if(length(clus) > 0){
  for(j in 1:length(clus)){
    pos <- which(mat3$clus %in% clus[j])
    mat4 <- mat3[pos,]
    chrinfo <- data.frame(table(mat4$chr))
    
    for(m in 1:nrow(chrinfo)){
      if(chrinfo$Freq[m] == 1){
        pos <- which(mat4$chr %in% chrinfo$Var1[m])
        mat4$clus[pos] <- "Other"
      }else{
        pos <- which(mat4$chr %in% chrinfo$Var1[m])
        mat5 <- mat4[pos,]
        mat5$st <- gsub(".*:","",mat5$ecDNA)
        mat5$en <- gsub(".*_","",mat5$st)
        mat5$st <- gsub("_.*","",mat5$st)
        mat5 <- mat5[order(as.numeric(as.character(mat5$st))),]
        mat5$type <- 0
        for(n in 2:nrow(mat5)){
          if((as.numeric(as.character(mat5$st[n])) - as.numeric(as.character(mat5$en[(n-1)]))) > 10000000){
            mat5$type[n] <- n
          }else{
            mat5$type[n] <- mat5$type[(n-1)]
          }
        }
        
        sui <- data.frame(table(mat5$type))
        for(n in 1:nrow(sui)){
          if(sui$Freq[n] == 1){
            pos <- which(mat5$type %in% sui$Var1[n])
            mat6 <- mat5[pos,]
            pos <- which(mat4$ecDNA %in% mat6$ecDNA)
            mat4$clus[pos] <- "Other"
          }else{
            pos <- which(mat5$type %in% sui$Var1[n])
            mat6 <- mat5[pos,]
            pos <- which(mat4$ecDNA %in% mat6$ecDNA)
            mat4$clus[pos] <- paste0(mat4$clus[pos],"_",chrinfo$Var1[m],"_",sui$Var1[n])
          }
        }
      }
    }
    mat.new <- rbind(mat.new,mat4)
    
  }
  }
}

mat.new1 <- NULL
for(i in 1:length(cellname)){
  pos <- which(mat.new$cellname %in% cellname[i])
  mattmp <- mat.new[pos,]
  mattmp$newclus <- "Other"
  pos <- which(mattmp$clus %in% c("Other","other"))
  mattmp1 <- mattmp[pos,]
  mattmp2 <- mattmp[-pos,]
  if(nrow(mattmp2) > 0){
    clus <- unique(mattmp2$clus)
    for(j in 1:length(clus)){
      pos <- which(mattmp2$clus %in% clus[j])
      mattmp2$newclus[pos] <- j
    }
  }
  mattmp3 <- rbind(mattmp1,mattmp2)
  mat.new1 <- rbind(mat.new1,mattmp3)
}

mat7 <- mat.new1[,c(1,2,6)]
colnames(mat7)[3] <- "clus"
siz <- NULL
for(i in 1:length(cellname)){
  pos <- which(mat7$cellname %in% cellname[i])
  mat8 <- mat7[pos,]
  pos1 <- which(mat8$clus %in% c("other","Other"))
  mat9 <- mat8[-pos,]
  siz <- c(siz,length(pos1) + length(unique(mat9$clus)))
}

siz1 <- c(764,254,118,172,1061,128,162,6,17,87,1029,1140,39,236,175,218,23,900,285,226,277,17,1,211,155,28,203,26,406,225,16)
summary(siz1)
summary(siz)
tes <- data.frame(cond = c(rep("Before",31),rep("After",31)),value = c(siz1,siz))
library(ggplot2)
ggplot(tes,aes(x= factor(cond,levels = c("Before","After")),y = value))+geom_boxplot()+theme_classic()+labs(x="",y="Number")

mat7$all <- paste0(mat7$cellname,":",mat7$clus)
mat7 <- mat7[order(mat7$cellname,mat7$all),]
write.table(mat7,file = "/Users/shangliu/01.terms/03.ecDNA/NC/ecDNA.set.all.new.xls",sep = "\t",quote=F,row.names = F)
pos <- which(mat2$cellname %in% cellname[8])
#matt <- mat2[pos,]

##########
mat7$all <- paste0(mat7$cellname,":",mat7$ecDNA)
pos <- which(mat7$clus %in% c("other","Other"))
matt <- mat7
alls <- unique(as.character(matt$all))
mit <- matrix(nrow = length(alls),ncol = 2)
for(i in 1:length(alls)){
  pos <- which(matt$all %in% alls[i])
  mat3 <- matt[pos,]
  siz <- 0
  for(j in 1:nrow(mat3)){
    tmp1 <- strsplit(mat3$ecDNA[j],split = ":")[[1]]
    tmp2 <- strsplit(tmp1[2],split = "_")[[1]]
    tmp3 <- as.numeric(as.character(tmp2[2])) - as.numeric(as.character(tmp2[1]))
    siz <- siz + tmp3
  }
  mit[i,1] <- alls[i]
  mit[i,2] <- siz
}
tmpd1 <- as.numeric(as.character(mit[,2])) / 5000000

mat7$all <- paste0(mat7$cellname,":",mat7$clus)
pos <- which(mat7$clus %in% c("other","Other"))
matt <- mat7[-pos,]
alls <- unique(as.character(matt$all))
mit <- matrix(nrow = length(alls),ncol = 2)
for(i in 1:length(alls)){
  pos <- which(matt$all %in% alls[i])
  mat3 <- matt[pos,]
  siz <- 0
  for(j in 1:nrow(mat3)){
    tmp1 <- strsplit(mat3$ecDNA[j],split = ":")[[1]]
    tmp2 <- strsplit(tmp1[2],split = "_")[[1]]
    tmp3 <- as.numeric(as.character(tmp2[2])) - as.numeric(as.character(tmp2[1]))
    siz <- siz + tmp3
  }
  mit[i,1] <- alls[i]
  mit[i,2] <- siz
}
tmpd2 <- as.numeric(as.character(mit[,2])) / 5000000

tmpd <- c(tmpd1,tmpd2)
summary(tmpd)
plot(density(c(tmpd,30)),main = "Length of ecDNA (M)")




