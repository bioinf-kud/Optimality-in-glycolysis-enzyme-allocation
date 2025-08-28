dir<-"-path to folder-/figures_and_source_data"
head1<-"/figureS7/data/"
head<-"/figureS7/data/group_gene_list/TCGA-"
tail<-"_gene_list_Tumor.tsv"
kegg_genes<-read.csv(paste0(dir,head1,"kegg_pathway_genes.csv"),header=TRUE)
gene_names<-read.table(paste0(dir,head1,"TCGA_name_to_ENSG.txt"),header=FALSE,sep='\t')
anno<-read.table(paste0(dir,head1,"ENSG_to_KEGG_HSA.txt"),sep='\t',header=TRUE)
data1<-read.table(paste0(dir,head,"BRCA",tail),header=TRUE,sep='\t')
data2<-read.table(paste0(dir,head,"COAD",tail),header=TRUE,sep='\t')
data3<-read.table(paste0(dir,head,"HNSC",tail),header=TRUE,sep='\t')
data4<-read.table(paste0(dir,head,"KIRC",tail),header=TRUE,sep='\t')
data5<-read.table(paste0(dir,head,"LUAD",tail),header=TRUE,sep='\t')
data6<-read.table(paste0(dir,head,"LUSC",tail),header=TRUE,sep='\t')
data7<-read.table(paste0(dir,head,"PRAD",tail),header=TRUE,sep='\t')
data8<-read.table(paste0(dir,head,"THCA",tail),header=TRUE,sep='\t')

names(anno)<-c("Column","Gene_ID")
data1<-merge(data1,anno,by="Column")
data2<-merge(data2,anno,by="Column")
data3<-merge(data3,anno,by="Column")
data4<-merge(data4,anno,by="Column")
data5<-merge(data5,anno,by="Column")
data6<-merge(data6,anno,by="Column")
data7<-merge(data7,anno,by="Column")
data8<-merge(data8,anno,by="Column")
data1<-subset(data1,select=-Column)
data2<-subset(data2,select=-Column)
data3<-subset(data3,select=-Column)
data4<-subset(data4,select=-Column)
data5<-subset(data5,select=-Column)
data6<-subset(data6,select=-Column)
data7<-subset(data7,select=-Column)
data8<-subset(data8,select=-Column)
anno<-merge(anno,gene_names,by.x="Column",by.y="V2")
names(anno)<-c("ENSEMBL","Gene_ID","Gene_name")
data1 <- aggregate(. ~ Gene_ID, data = data1, FUN = function(x) {
  sum(x) 
})
data2 <- aggregate(. ~ Gene_ID, data = data2, FUN = function(x) {
  sum(x) 
})
data3 <- aggregate(. ~ Gene_ID, data = data3, FUN = function(x) {
  sum(x) 
})
data4 <- aggregate(. ~ Gene_ID, data = data4, FUN = function(x) {
  sum(x) 
})
data5 <- aggregate(. ~ Gene_ID, data = data5, FUN = function(x) {
  sum(x) 
})
data6 <- aggregate(. ~ Gene_ID, data = data6, FUN = function(x) {
  sum(x) 
})
data7 <- aggregate(. ~ Gene_ID, data = data7, FUN = function(x) {
  sum(x) 
})
data8 <- aggregate(. ~ Gene_ID, data = data8, FUN = function(x) {
  sum(x) 
})
row.names(data1)<-data1$Gene_ID
row.names(data2)<-data2$Gene_ID
row.names(data3)<-data3$Gene_ID
row.names(data4)<-data4$Gene_ID
row.names(data5)<-data5$Gene_ID
row.names(data6)<-data6$Gene_ID
row.names(data7)<-data7$Gene_ID
row.names(data8)<-data8$Gene_ID
data1$samp="BRCA"
data2$samp="COAD"
data3$samp="HNSC"
data4$samp="KIRC"
data5$samp="LUAD"
data6$samp="LUSC"
data7$samp="PRAD"
data8$samp="THCA"
matrix<-rbind(data1,data2)
matrix<-rbind(matrix,data3)
matrix<-rbind(matrix,data4)
matrix<-rbind(matrix,data5)
matrix<-rbind(matrix,data6)
matrix<-rbind(matrix,data7)
matrix<-rbind(matrix,data8)
matrix<-merge(matrix,anno,by="Gene_ID")
all<-matrix
matrix<-matrix[matrix$Gene_ID %in% kegg_genes[kegg_genes$Pathway_ID=="hsa03010",]$Gene_ID,]
ggplot(matrix, aes(x=samp, y = Spearman_Correlation)) +
  geom_violin(
    fill = '#AFCBE3', 
    color = "#4472C4", 
    alpha = 0.6, 
    draw_quantiles = c(0.5)  
  ) +
  scale_color_manual(values = DEcolor) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#4472C4", size = 1)+
  ggtitle("Correlation between expression of ribosome genes \nand EOI of glycolysis in human tumor samples") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),  
    axis.text.y = element_text(size = 20),  
    axis.title = element_text(size = 20),  
    plot.title = element_text(size = 22, hjust = 0.5, vjust = 1), 
    plot.title.position = "plot",  
    panel.grid = element_blank(),  
    panel.background = element_rect(fill = "white", color = NA),  
    panel.border = element_rect(color = "black", fill = NA, size = 1)  
  ) +
  labs(y = expression("Spearman's " ~ rho)) +
  xlab("") +
  ylim(-1, 1)
matrix<-all
matrix<-matrix[matrix$Gene_ID %in% kegg_genes[kegg_genes$Pathway_ID=="hsa00190",]$Gene_ID,]
ggplot(matrix, aes(x=samp, y = Spearman_Correlation)) +
  geom_violin(
    fill = '#AFCBE3', 
    color = "#4472C4", 
    alpha = 0.6, 
    draw_quantiles = c(0.5)  
  ) +
  scale_color_manual(values = DEcolor) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#4472C4", size = 1)+
  ggtitle("Correlation between expression of OXPHOS genes \nand EOI of glycolysis in human tumor samples") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),  
    axis.text.y = element_text(size = 20),  
    axis.title = element_text(size = 20),  
    plot.title = element_text(size = 22, hjust = 0.5, vjust = 1), 
    plot.title.position = "plot",  
    panel.grid = element_blank(),  
    panel.background = element_rect(fill = "white", color = NA),  
    panel.border = element_rect(color = "black", fill = NA, size = 1)  
  ) +
  labs(y = expression("Spearman's " ~ rho)) +
  xlab("") +
  ylim(-1, 1)






