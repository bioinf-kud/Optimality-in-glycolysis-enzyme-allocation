library(dplyr)
library(ggplot2)
library(ggrepel) 

dir<-"-path to folder-/figures_and_source_data"
file <- paste0(dir, "/figure 5/data/corr_P.tsv")
data <- read.table(file, header = TRUE, sep = "\t")

ggplot(data, aes(x = rho, y = -log10(p), color = species, label = type)) +
  geom_point(alpha = 0.5, size = 5) +
  geom_text_repel(
    data = data[data$p < 0.05, ],
    size = 5,
    max.overlaps = 30,        
    box.padding = 0.5,       
    segment.color = "grey50"  
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#824098", size = 1, alpha=0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#824098", size = 1,alpha=0.7) +
  scale_color_manual(
    values = c(
      "H.sapiens" = "#4472C4", 
      "E.coli" = "#ED7D31", 
      "S.cerevisiae" = "#70AD47"
    ), 
    labels = c(
      "H.sapiens" = expression(italic("H.sapiens")), 
      "E.coli" = expression(italic("E.coli")), 
      "S.cerevisiae" = expression(italic("S.cerevisiae"))
    ), 
    name = "species"
  ) +
  labs(
    x = expression("Spearman's " ~ rho),
    y = expression(-log[10]*"(p-value)"),
    color = "Experimental condition",
    title = "Correlation between glycolysis enzyme\n fraction and EOI of glycolysis"
  ) +
  theme(
    axis.text.x = element_text(hjust = 1, size = 20),  
    axis.text.y = element_text(size = 20),  
    axis.title = element_text(size = 20),  
    plot.title = element_text(size = 22, hjust = 0.5, vjust = 1), 
    plot.title.position = "plot",  
    panel.grid = element_blank(),  
    panel.background = element_rect(fill = "white", color = NA),  
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    legend.background = element_rect(fill = "transparent", color = NA), 
    legend.text = element_text(size = 16),        
    legend.title = element_text(size = 18), 
    legend.position = c(0.85,0.80)
  )

