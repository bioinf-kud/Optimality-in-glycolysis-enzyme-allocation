library(dplyr)
library(ggplot2)

dir<-"-path to folder-/figures_and_source_data"
file <- paste0(dir, "/figureS9/data/Ecoli_sum_corr.tsv")
aggrt <- read.table(file, header = TRUE, sep = "\t",)

corr_stats <- aggrt %>%
  filter(type != "Carbon source") %>%
  group_by(type) %>%
  summarise(
    rho = cor(Glyco_Genes_Sum, Corr, method = "spearman", use = "complete.obs"),
    p = cor.test(Glyco_Genes_Sum, Corr, method = "spearman")$p.value
  ) %>%
  mutate(
    Rho_label = paste0(type, " (", round(rho, 2), ", ", formatC(p, format="e", digits=2), ")")
  )

aggrt <- aggrt %>%
  left_join(corr_stats %>% select(type, Rho_label), by="type")

aggrt$Rho_label[aggrt$type == "Carbon source"] <- "Carbon source"

group_colors <- c(
  "C-limitation" = "#4472C4", 
  "A-limitation" = "#ED7D31", 
  "R-limitation" = "#70AD47", 
  "Carbon source" = "#824098"
)
new_group_colors <- setNames(
  c(group_colors[corr_stats$type], group_colors["Carbon source"]),
  c(corr_stats$Rho_label, "Carbon source")
)
ggplot(aggrt, aes(x = Glyco_Genes_Sum, y = Corr, color = Rho_label, label = Description)) +
  geom_point(alpha = 0.5,size=4) +
  geom_smooth(
    data = subset(aggrt, type != "Carbon source"), 
    method = "lm", se = FALSE
  ) +
  geom_text(
    data = subset(aggrt, type == "Carbon source"),
    size = 5, hjust = 0, vjust = 1, check_overlap = TRUE
  )+
  scale_color_manual(values = new_group_colors, name = "Group (R²)")+
  labs(
    x = "Mass fraction of glycolysis enzyme",
    y = "EOI of glycolysis",
    color = "Experimental condition",
    title = expression(
      atop(
        "Expression of glycolysis enzyme and",
        "EOI of glycolysis in "* italic("E.coli")*" samples"
      )
    )
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
    legend.position = c(0.72,0.85)
  )

