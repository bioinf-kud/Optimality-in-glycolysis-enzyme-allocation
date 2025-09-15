library(clusterProfiler)
library(enrichplot)
library(gridExtra)
library(grid)

dir<-"-path to folder-/figures_and_source_data"
head1<-"/figure 5/data/"
head<-"/figure 5/data/group_gene_list/TCGA-"
tail<-"_gene_list_Tumor.tsv"
pathway<-read.csv(paste0(dir,head1,"pathway_nodes_with_categories.csv"))
anno<-read.table(paste0(dir,head1,"ENSG_to_KEGG_HSA.txt"),sep='\t',header=TRUE)
data2<-read.table(paste0(dir,head,"COAD",tail),header=TRUE,sep='\t')
names(anno)<-c("Column","Gene_ID")
data2<-merge(data2,anno,by="Column")
data2<-subset(data2,select=-Column)
data2 <- aggregate(. ~ Gene_ID, data = data2, FUN = function(x) {
  sum(x) 
})
row.names(data2)<-data2$Gene_ID
data2<-subset(data2,select=-Gene_ID)








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
plot_top3_pathway_gsea(gseKEGG2, "COAD", pathway)
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
