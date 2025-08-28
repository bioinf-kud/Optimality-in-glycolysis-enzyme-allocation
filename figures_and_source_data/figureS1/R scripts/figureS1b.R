library(ggplot2)

dir<-"-path to folder-/figures_and_source_data"
head <- "/figureS1/data/avgs.tsv"
data <- read.table(paste0(dir, head), header = TRUE, sep = '\t')

spearman <- cor.test(log(data$tran_avg), log(data$prot_avg), method = "spearman")
rho <- round(spearman$estimate, 3)
pval <- signif(spearman$p.value, 2)

ggplot(data, aes(x = log(tran_avg), y = log(prot_avg))) +
  geom_point(alpha = 0.5, color = "#4472C4") +
  labs(
    title = "Correlation of mean protein and transcript\nabundance for human metabolic enzymes",
    x = "Log mean TPM, TCGA transcriptomics",
    y = "Log mean abundance, CCLE proteomics"
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
  annotate("text", 
           x = min(log(data$tran_avg), na.rm = TRUE), 
           y = max(log(data$prot_avg), na.rm = TRUE) + 0.5, 
           label = paste0("Spearman's \u03C1 = ", rho, "\nP < 1e-10 "), 
           hjust = -0.01, vjust = 1, size = 7, fontface = "plain", color = "black")
