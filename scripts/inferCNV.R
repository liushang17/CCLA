library(infercnv)
library(HiddenMarkov)
library(Seurat)
args <- commandArgs(T)

# input ref data 
data <- read.table("/mydata/RPE1_count_mtx.tsv.txt",sep="\t",header=T,row.names=1)
meta <- data.frame(Cluster=rep("RPE1",dim(data)[2]),row.names = colnames(data))
ob1 <- CreateSeuratObject(data,meta.data = meta)

data <- read.table("/mydata/HK2_count_mtx.tsv.txt",sep="\t",header=T,row.names=1)
meta <- data.frame(Cluster=rep("HK2",dim(data)[2]),row.names = colnames(data))
ob2 <- CreateSeuratObject(data,meta.data = meta)

ob1 <- merge(ob1,ob2)

# input query data and merge into ref data 
i<-as.numeric(args[1])
path <- "/mydata/RNA_counts_data"
path <- list.files(path,pattern=".rds",full.names=T)
ob <- readRDS(path[i])
pro <- gsub("_seurat.rds","",basename(path[i])) 
ob$Cluster <- ob$seurat_clusters
ob <-merge(ob,ob1)
data <- as.matrix(GetAssayData(ob,slot="counts"))
data <- data[which(rowSums(data)>0),]
meta <- data.frame(Sample = ob$Cluster,row.names = colnames(ob))

# creat outdir 
outdir<-"/jdfssz1/ST_SUPERCELLS/P21Z10200N0134/USER/zhangyuanhang1/CRC/CCLA/update_0425/01.result/umap.123/20230322/test/infercnv/"
dir.create(paste0(outdir,pro))
outdir<-paste0(paste0(outdir,pro))

# infercnv 
bin<-400
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=data,annotations_file=meta,
                                    delim="\t",gene_order_file="/jdfssz1/ST_SUPERCELLS/P21Z10200N0134/USER/zhangyuanhang1/CRC/03.CNV/gene.pos_4.txt",
                                    ref_group_names=c("RPE1","HK2"))

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir=outdir,
                             window_length = bin, #default 101
                             cluster_by_groups=T,analysis_mode = "subclusters",
                             denoise=T,
                             HMM=T)
