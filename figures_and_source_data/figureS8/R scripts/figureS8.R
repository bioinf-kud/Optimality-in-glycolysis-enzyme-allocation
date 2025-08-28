library(clusterProfiler)
library(enrichplot)
library(gridExtra)
library(grid)


dir<-"-path to folder-/figures_and_source_data"
head1<-"/figureS8/data/"
head<-"/figureS8/data/group_gene_list/TCGA-"
tail<-"_gene_list_Tumor.tsv"
pathway<-read.csv(paste0(dir,head1,"pathway_nodes_with_categories.csv"))
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
data1<-subset(data1,select=-Gene_ID)
data2<-subset(data2,select=-Gene_ID)
data3<-subset(data3,select=-Gene_ID)
data4<-subset(data4,select=-Gene_ID)
data5<-subset(data5,select=-Gene_ID)
data6<-subset(data6,select=-Gene_ID)
data7<-subset(data7,select=-Gene_ID)
data8<-subset(data8,select=-Gene_ID)







data<-data1
genelist<-data[,1]
names(genelist) = row.names(data)
genelist = sort(genelist, decreasing = TRUE)
gseKEGG1 <- gseKEGG(geneList     = genelist,
                   organism     = 'hsa',
                   minGSSize    = 100,
                   maxGSSize    = 500,
                   pvalueCutoff = 0.05,
                   verbose      = FALSE)
gseKEGG<-as.data.frame(gseKEGG1)
gseKEGG<-subset(gseKEGG,select=c("ID","enrichmentScore"))
names(gseKEGG)<-c("ID","BRCA")
allgseKEGG<-gseKEGG

data<-data2
genelist<-data[,1]
names(genelist) = row.names(data)
genelist = sort(genelist, decreasing = TRUE)
gseKEGG2 <- gseKEGG(geneList     = genelist,
                   organism     = 'hsa',
                   minGSSize    = 100,
                   maxGSSize    = 500,
                   pvalueCutoff = 0.05,
                   verbose      = FALSE)
gseKEGG<-as.data.frame(gseKEGG2)
gseKEGG<-subset(gseKEGG,select=c("ID","enrichmentScore"))
names(gseKEGG)<-c("ID","COAD")
allgseKEGG<-merge(allgseKEGG,gseKEGG,by="ID",all = TRUE)

data<-data3
genelist<-data[,1]
names(genelist) = row.names(data)
genelist = sort(genelist, decreasing = TRUE)
gseKEGG3 <- gseKEGG(geneList     = genelist,
                   organism     = 'hsa',
                   minGSSize    = 100,
                   maxGSSize    = 500,
                   pvalueCutoff = 0.05,
                   verbose      = FALSE)
gseKEGG<-as.data.frame(gseKEGG3)
gseKEGG<-subset(gseKEGG,select=c("ID","enrichmentScore"))
names(gseKEGG)<-c("ID","HNSC")
allgseKEGG<-merge(allgseKEGG,gseKEGG,by="ID",all = TRUE)

data<-data4
genelist<-data[,1]
names(genelist) = row.names(data)
genelist = sort(genelist, decreasing = TRUE)
gseKEGG4 <- gseKEGG(geneList     = genelist,
                   organism     = 'hsa',
                   minGSSize    = 100,
                   maxGSSize    = 500,
                   pvalueCutoff = 0.05,
                   verbose      = FALSE)
gseKEGG<-as.data.frame(gseKEGG4)
gseKEGG<-subset(gseKEGG,select=c("ID","enrichmentScore"))
names(gseKEGG)<-c("ID","KIRC")
allgseKEGG<-merge(allgseKEGG,gseKEGG,by="ID",all = TRUE)

data<-data5
genelist<-data[,1]
names(genelist) = row.names(data)
genelist = sort(genelist, decreasing = TRUE)
gseKEGG5 <- gseKEGG(geneList     = genelist,
                   organism     = 'hsa',
                   minGSSize    = 100,
                   maxGSSize    = 500,
                   pvalueCutoff = 0.05,
                   verbose      = FALSE)
gseKEGG<-as.data.frame(gseKEGG5)
gseKEGG<-subset(gseKEGG,select=c("ID","enrichmentScore"))
names(gseKEGG)<-c("ID","LUAD")
allgseKEGG<-merge(allgseKEGG,gseKEGG,by="ID",all = TRUE)

data<-data6
genelist<-data[,1]
names(genelist) = row.names(data)
genelist = sort(genelist, decreasing = TRUE)
gseKEGG6 <- gseKEGG(geneList     = genelist,
                   organism     = 'hsa',
                   minGSSize    = 100,
                   maxGSSize    = 500,
                   pvalueCutoff = 0.05,
                   verbose      = FALSE)
gseKEGG<-as.data.frame(gseKEGG6)
gseKEGG<-subset(gseKEGG,select=c("ID","enrichmentScore"))
names(gseKEGG)<-c("ID","LUSC")
allgseKEGG<-merge(allgseKEGG,gseKEGG,by="ID",all = TRUE)

data<-data7
genelist<-data[,1]
names(genelist) = row.names(data)
genelist = sort(genelist, decreasing = TRUE)
gseKEGG7 <- gseKEGG(geneList     = genelist,
                   organism     = 'hsa',
                   minGSSize    = 100,
                   maxGSSize    = 500,
                   pvalueCutoff = 0.05,
                   verbose      = FALSE)
gseKEGG<-as.data.frame(gseKEGG7)
gseKEGG<-subset(gseKEGG,select=c("ID","enrichmentScore"))
names(gseKEGG)<-c("ID","PRAD")
allgseKEGG<-merge(allgseKEGG,gseKEGG,by="ID",all = TRUE)

data<-data8
genelist<-data[,1]
names(genelist) = row.names(data)
genelist = sort(genelist, decreasing = TRUE)
gseKEGG8 <- gseKEGG(geneList     = genelist,
                   organism     = 'hsa',
                   minGSSize    = 100,
                   maxGSSize    = 500,
                   pvalueCutoff = 0.05,
                   verbose      = FALSE)
gseKEGG<-as.data.frame(gseKEGG8)
gseKEGG<-subset(gseKEGG,select=c("ID","enrichmentScore"))
names(gseKEGG)<-c("ID","THCA")
allgseKEGG<-merge(allgseKEGG,gseKEGG,by="ID",all = TRUE)
names(pathway)<-c("ID","name")


plot_top3_pathway_gsea <- function(gsea_result, cancer_name, pathway_df, colors = c("#4472C4", "#ED7D31", "#70AD47")) {
  df <- as.data.frame(gsea_result)
  df_sub <- df[df$ID %in% pathway_df$ID, ]
  top3 <- head(df_sub[order(-df_sub$NES), ], 3)
  top3_ids <- top3$ID
  top3_desc <- top3$Description
  if(length(top3_ids) == 0){
    message("No valid pathway found for ", cancer_name)
    return(NULL)
  }
  title_str <- paste0("Top pathways enriched in genes positively\ncorrelated with EOI of glycolysis in TCGA-", cancer_name)
  gseaplot2(
    gsea_result,
    geneSetID = top3_ids,
    color = colors[seq_along(top3_ids)],
    pvalue_table = TRUE,
    title = title_str,
    base_size = 14,
    ES_geom = "line",
  )
}


plot_top3_pathway_gsea(gseKEGG1, "BRCA", pathway)
grid::grid.rect(
  x = 0.02, y = 0.13,    
  width = 0.04, height = 0.5,  
  gp = grid::gpar(fill = "white", col = NA), 
  just = "center"
)
grid::grid.text(
  "Correlation with EOI of glycolysis",
  x = 0.03, y = 0.22,
  rot = 90,
  gp = grid::gpar(fontsize = 12, col = "black")
)

plot_top3_pathway_gsea(gseKEGG3, "HNSC", pathway)
grid::grid.rect(
  x = 0.02, y = 0.13,    
  width = 0.04, height = 0.5,  
  gp = grid::gpar(fill = "white", col = NA), 
  just = "center"
)
grid::grid.text(
  "Correlation with EOI of glycolysis",
  x = 0.03, y = 0.22,
  rot = 90,
  gp = grid::gpar(fontsize = 12, col = "black")
)
plot_top3_pathway_gsea(gseKEGG4, "KIRC", pathway)
grid::grid.rect(
  x = 0.02, y = 0.13,    
  width = 0.04, height = 0.5,  
  gp = grid::gpar(fill = "white", col = NA), 
  just = "center"
)
grid::grid.text(
  "Correlation with EOI of glycolysis",
  x = 0.03, y = 0.22,
  rot = 90,
  gp = grid::gpar(fontsize = 12, col = "black")
)
plot_top3_pathway_gsea(gseKEGG5, "LUAD", pathway)
grid::grid.rect(
  x = 0.02, y = 0.13,    
  width = 0.04, height = 0.5,  
  gp = grid::gpar(fill = "white", col = NA), 
  just = "center"
)
grid::grid.text(
  "Correlation with EOI of glycolysis",
  x = 0.03, y = 0.22,
  rot = 90,
  gp = grid::gpar(fontsize = 12, col = "black")
)
plot_top3_pathway_gsea(gseKEGG6, "LUSC", pathway)
grid::grid.rect(
  x = 0.02, y = 0.13,    
  width = 0.04, height = 0.5,  
  gp = grid::gpar(fill = "white", col = NA), 
  just = "center"
)
grid::grid.text(
  "Correlation with EOI of glycolysis",
  x = 0.03, y = 0.22,
  rot = 90,
  gp = grid::gpar(fontsize = 12, col = "black")
)
plot_top3_pathway_gsea(gseKEGG7, "PRAD", pathway)
grid::grid.rect(
  x = 0.02, y = 0.13,    
  width = 0.04, height = 0.5,  
  gp = grid::gpar(fill = "white", col = NA), 
  just = "center"
)
grid::grid.text(
  "Correlation with EOI of glycolysis",
  x = 0.03, y = 0.22,
  rot = 90,
  gp = grid::gpar(fontsize = 12, col = "black")
)

