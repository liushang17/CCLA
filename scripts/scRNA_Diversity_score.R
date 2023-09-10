library(Seurat)
library(ggplot2)
library(dplyr)
library(forcats)
library(ggpubr)

load("mydata/allcellines_scRNA.Rdata")

#compute diversity score
pca <- ob@reductions$pca@cell.embeddings
pca <- pca[,1:30]
ITH <- c()
for (i in unique(ob$Alias_update)) {
  temp_mat <- pca[which(ob$Alias_update==i),]
  cutoff1 <- apply(temp_mat, 2, function(x) mean(x)-3*sd(x))
  cutoff2 <- apply(temp_mat, 2, function(x) mean(x)+3*sd(x))
  cell_bool <- apply(temp_mat, 1, function(x) sum(x[1:3] < cutoff1[1:3]) < 3 & sum(x[1:3] > cutoff2[1:3]) < 3)
  temp_mat <- temp_mat[cell_bool,]
  mean_m<-apply(temp_mat, 2, function(y) mean(y))
  diversity <- mean(apply(temp_mat, 1, function(x) sqrt(sum((x-mean_m)^2))))
  ITH <- c(ITH,diversity)
}
names(ITH) <-unique(ob$Alias_update)
ITH<-as.data.frame(ITH)
colnames(ITH)<-"Diversity_score"
ITH$Alias_update <-rownames(ITH)

# input Continuous_Discrete.csv
dc <- read.csv("Continuous_Discrete.csv")
data_plot <- left_join(ITH,dc)
#write.csv(data_plot,"allcellines_ITH_DC.csv")
data_plot  <- data_plot[-which(data_plot$Cell_line %in% c("Caco2_2","MDAMB231_2","SCC4_2")),]
# vlnplot  
p1<- ggplot(data_plot,aes(Cluster_status,Diversity_score,fill=Cluster_status))+
  geom_violin(trim = FALSE)+
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
  geom_signif(comparisons = list(c("Continuous", "Discrete")),
              map_signif_level=F,
              textsize=6,test="wilcox.test",test.args = "less",step_increase=0.2)
#barplot 
p2<- data_plot %>%
  mutate(Alias_update = fct_reorder(Alias_update,Diversity_score)) %>%
  ggplot(aes(Alias_update,Diversity_score,fill=Cluster_status)) +
  geom_bar(stat = "identity")+
  xlab("") +
  theme_bw()+
  theme(panel.grid=element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = 6),
        #axis.title.y=element_text(size=8,family = "Arial")
        #axis.line = element_line(size = 0.8, colour = "black")
  )+
  #scale_y_continuous(limits = c(10,30),breaks=seq(10,30,10),expand = c(0,0))+
  scale_fill_manual(values = c("#d0d1e6","#0570b0"))

pdf("Diversity_VlnPlot.pdf",4,4)
p1
dev.off()
pdf("Diversity_barplot.pdf",7,4)
p2
dev.off()




