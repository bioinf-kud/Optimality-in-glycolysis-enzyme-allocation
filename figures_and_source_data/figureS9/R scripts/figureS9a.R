library(dplyr)
library(ggplot2)

dir<-"-path to folder-/figures_and_source_data"
file <- paste0(dir, "/figureS9/data/CCLE_sum_corr.tsv")
aggrt <- read.table(file, header = TRUE, sep = "\t",)
corr_stats <- aggrt %>%
  group_by(type) %>%
  summarise(
    rho = cor(Glyco_Genes_Sum, Corr, method = "spearman", use = "complete.obs"),
    p = cor.test(Glyco_Genes_Sum, Corr, method = "spearman")$p.value
  ) %>%
  mutate(
    Rho_label = paste0(type, " (", round(rho, 2), ", ", formatC(p, format="e", digits=2), ")")
  )
aggrt <- aggrt %>% left_join(corr_stats %>% select(type, Rho_label), by="type")
group_colors <- c(
  "LUNG" = "#4472C4", 
  "HAEMATOPOIETIC" = "#ED7D31", 
  "LARGE_INTESTINE" = "#70AD47", 
  "SKIN" = "#824098",
  "BREAST" = "#C77EB5",
  "PANCREAS" = "#FFD166",
  "OVARY" = "#B89C6F"
)
new_group_colors <- setNames(group_colors[corr_stats$type], corr_stats$Rho_label)

ggplot(aggrt, aes(x = Glyco_Genes_Sum, y = Corr, color = Rho_label)) +
  geom_point(alpha = 0.5, size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = new_group_colors, name = expression("Cell line source (" ~ rho ~", p-value)")) +
  facet_wrap(~ type, scales = "free") +
  labs(
    x = "Mass fraction of lycolysis enzyme",
    y = "EOI of glycolysis",
    color = "Cell line source",
    title = "Expression of glycolysis enzyme and EOI of glycolysis in human cell lines"
  ) +
  theme(
    axis.title.x = element_text(size = 16),  
    axis.text.y = element_text(size = 16),  
    axis.text.x = element_text(size = 16), 
    axis.title = element_text(size = 16),  
    plot.title = element_text(size = 22, hjust = 0.5, vjust = 1), 
    plot.title.position = "plot",  
    panel.grid = element_blank(),  
    panel.background = element_rect(fill = "white", color = NA),  
    panel.border = element_rect(color = "black", fill = NA, size = 1),  
    legend.position = c(0.64, 0.17), 
    legend.background = element_rect(fill = "transparent", color = NA), 
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 20),
    strip.text = element_text(size = 16)   
  ) 