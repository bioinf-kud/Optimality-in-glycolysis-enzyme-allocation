library(ggplot2)
library(ggrepel)
dir<-"-path to folder-/figures_and_source_data"
head<-"/figure 1/data/"
tail<-".tsv"
CCLE<-read.table(paste0(dir,head,"CCLE",tail),header=TRUE,sep='\t')


p<-ggplot(CCLE) +
  geom_histogram(
    aes(x = SpearmanCorrelation), 
    binwidth = 0.010, 
    position = "identity", 
    alpha = 0.5, 
    fill = "#4472C4",      
    colour = "#4472C4"   
  ) +
  labs(
    title = expression(
      atop(
        "Correlation of metabolic enzyme abundance ",
        "between 378 cancer cell lines"
      )
    ),
    x = expression("Spearman's " ~ rho),
    y = "Number of cell line pairs"
  ) +
  xlim(0.6, 1) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 22, hjust = 0.5, vjust = 1),
    plot.title.position = "plot",
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )
ggsave(paste0(dir,"/figure 1/figures/","figure 1c.pdf"), p, width = 650/96, height = 700/96, units = "in", device = "pdf")

