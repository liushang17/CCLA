library(ggplot2)

indir <- "/Users/shangliu/01.terms/03.ecDNA/NC/08.Siggenes/"
files <- list.files(indir)
pos <- grep("_ecDNA_in_markers.txt",files)
files1 <- files[pos]
pos <- grep("_ecDNA_in_tfs.txt",files)
files2 <- files[pos]


mat <- NULL
for(i in 1:length(files1)){
  infile <- paste0(indir,files1[i])
  mattmp <- read.table(infile,sep = "\t",header = T)
  if(nrow(mattmp) > 0){
  mattmp$tmp <- gsub("_ecDNA_in_markers.txt","",files1[i])
  mat <- rbind(mat,mattmp)
  }
}

tfs <- NULL
for(i in 1:length(files2)){
  infile <- paste0(indir,files2[i])
  mattmp <- read.table(infile,sep = "\t")
  if(nrow(mattmp) > 0){
    mattmp$tmp <- gsub("_ecDNA_in_tfs.txt","",files2[i])
    tfs <- rbind(tfs,mattmp)
  }
}

mat$all <- paste0(mat$tmp,"_",mat$gene)
tfs$all <- paste0(tfs$tmp,"_",tfs$V1)

pos <- which(mat$all %in% tfs$all)
mat1 <- mat[pos,]

mit <- data.frame(table(tfs$tmp))
mit$type <- "TFs"
pos <- which(mit$Freq > 30)
mit$Freq[pos] <- 30
ggplot(mit,aes(x= type,y = Freq))+geom_boxplot()+theme_classic()

genes <- read.table("/Users/shangliu/01.terms/03.ecDNA/NC/05.getgene/SNB75_CL100169085_L02_C6_gene.txt",sep = "\t")
mit <- read.table("/Users/shangliu/01.terms/03.ecDNA/NC/08.ecDNA_cell_mtx/SNB75_CL100169085_L02_C6.ecDNA_cell_mtx.txt",sep = "\t",header = F)

pos <- which(genes$V2 %in% "SOX10")
genes[pos,]


