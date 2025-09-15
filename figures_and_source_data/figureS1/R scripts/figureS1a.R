library(ggplot2)
dir<-"-path to folder-/figures_and_source_data"
head<-"/figureS1/data/"
tail<-".tsv"
CCLE<-read.table(paste0(dir,head,"pro_trans",tail),header=TRUE,sep='\t')


ggplot(CCLE) +
  geom_histogram(
    aes(x = SpearmanCorrelation), 
    binwidth = 0.020, 
    position = "identity", 
    alpha = 0.5, 
    fill = "#4472C4",      
    colour = "#4472C4"   
  ) +
  labs(
    title = expression(
      atop(
        "Correlation of metabolic enzyme abundance",
        "between proteomes and transcriptomes"
      )
    ),
    x = expression("Spearman's " ~ rho),
    y = "Number of proteome-transcriptom pairs"
  ) +
  xlim(0, 1) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 22, hjust = 0.5, vjust = 1),
    plot.title.position = "plot",
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )+
  scale_y_continuous(
  labels = function(x) parse(text = paste0(x/1e5, " %*% 10^5"))
)
