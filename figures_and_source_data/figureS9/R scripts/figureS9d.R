library(dplyr)
library(ggplot2)

dir<-"-path to folder-/figures_and_source_data"
file <- paste0(dir, "/figureS9/data/yeast_sum_corr.tsv")
aggrt <- read.table(file, header = TRUE, sep = "\t",)

spearman <- cor.test(aggrt$molecule.fraction.of.glycolysis.enzyme, aggrt$correlation, method = "spearman")
rho <- round(spearman$estimate, 3)
pval <- signif(spearman$p.value, 2)

ggplot(aggrt, aes(x = molecule.fraction.of.glycolysis.enzyme, y = correlation)) +
  geom_point(alpha = 0.5, color = "#4472C4",size=5) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    x = "Molecular fraction of glycolysis enzyme ",
    y = "EOI of glycolysis",
    color = "Experimental condition",
    title = expression(
      atop(
        "Expression of glycolysis enzyme and",
        "EOI of glycolysis in "* italic("S.cerevisiae")*" samples"
      )
    )
  ) +
  theme(
    axis.text.x = element_text(hjust = 1, size = 20),  
    axis.text.y = element_text(size = 20),  
    axis.title = element_text(size = 20, face = "plain"),  
    plot.title = element_text(size = 22, hjust = 0.5, vjust = 1, face = "plain"), 
    plot.title.position = "plot",  
    panel.grid = element_blank(),  
    panel.background = element_rect(fill = "white", color = NA),  
    panel.border = element_rect(color = "black", fill = NA, size = 1)  
  ) +
  annotate(
    "text", 
    x = max(aggrt$molecule.fraction.of.glycolysis.enzyme)/2, 
    y = max(aggrt$correlation), 
    label = paste0("Spearman's \u03C1 = ", rho, "\nP<1e-10 "), 
    hjust = 0, vjust = 1 , size = 7, fontface = "plain", color = "black"
  )



