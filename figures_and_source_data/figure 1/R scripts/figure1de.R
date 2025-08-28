library(ggplot2)
library(ggrepel)
dir<-"-path to folder-/figures_and_source_data"
head<-"/figure 1/data/"
tail<-"_tumor.tsv"
BRCA<-read.table(paste0(dir,head,"BRCA",tail),header=TRUE,sep='\t')
COAD<-read.table(paste0(dir,head,"COAD",tail),header=TRUE,sep='\t')
HNSC<-read.table(paste0(dir,head,"HNSC",tail),header=TRUE,sep='\t')
KIRC<-read.table(paste0(dir,head,"KIRC",tail),header=TRUE,sep='\t')
LUAD<-read.table(paste0(dir,head,"LUAD",tail),header=TRUE,sep='\t')
LUSC<-read.table(paste0(dir,head,"LUSC",tail),header=TRUE,sep='\t')
PRAD<-read.table(paste0(dir,head,"PRAD",tail),header=TRUE,sep='\t')
THCA<-read.table(paste0(dir,head,"THCA",tail),header=TRUE,sep='\t')
merged<-rbind(BRCA,COAD)
merged<-rbind(merged,HNSC)
merged<-rbind(merged,KIRC)
merged<-rbind(merged,LUAD)
merged<-rbind(merged,LUSC)
merged<-rbind(merged,PRAD)
merged<-rbind(merged,THCA)

ggplot(merged, aes(x = Group, y = SpearmanCorrelation)) +
  geom_violin(
    fill = '#AFCBE3', 
    color = "#4472C4", 
    alpha = 0.6, 
    draw_quantiles = c(0.5)  
  ) +
  scale_color_manual(values = DEcolor) +
  ggtitle("Correlation of metabolic gene expression\nbetween human tumor samples") +
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
  ylim(0, 1)
  
  


dir<-"/Users/sunkai/Desktop/metabolomic_proj/文章图片及数据/figure 1"
head<-"/data/"
tail<-"_normal.tsv"
BRCA<-read.table(paste0(dir,head,"BRCA",tail),header=TRUE,sep='\t')
COAD<-read.table(paste0(dir,head,"COAD",tail),header=TRUE,sep='\t')
HNSC<-read.table(paste0(dir,head,"HNSC",tail),header=TRUE,sep='\t')
KIRC<-read.table(paste0(dir,head,"KIRC",tail),header=TRUE,sep='\t')
LUAD<-read.table(paste0(dir,head,"LUAD",tail),header=TRUE,sep='\t')
LUSC<-read.table(paste0(dir,head,"LUSC",tail),header=TRUE,sep='\t')
PRAD<-read.table(paste0(dir,head,"PRAD",tail),header=TRUE,sep='\t')
THCA<-read.table(paste0(dir,head,"THCA",tail),header=TRUE,sep='\t')
merged<-rbind(BRCA,COAD)
merged<-rbind(merged,HNSC)
merged<-rbind(merged,KIRC)
merged<-rbind(merged,LUAD)
merged<-rbind(merged,LUSC)
merged<-rbind(merged,PRAD)
merged<-rbind(merged,THCA)

ggplot(merged, aes(x = Group, y = SpearmanCorrelation)) +
  geom_violin(
    fill = '#AFCBE3', 
    color = "#4472C4", 
    alpha = 0.6, 
    draw_quantiles = c(0.5)  
  ) +
  scale_color_manual(values = DEcolor) +
  ggtitle("Correlation of metabolic gene expression\nbetween normal human tissues") +
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
  ylim(0, 1)

