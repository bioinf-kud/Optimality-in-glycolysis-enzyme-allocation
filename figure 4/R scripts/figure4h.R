library(ggplot2)
library(ggrepel)
library(ggbeeswarm)

dir<-"-path to folder-/figures_and_source_data"
head <- "/figure 4/data/corr"
tail <- ".csv"
df <- read.csv(paste0(dir, head, tail), header = TRUE)


species_colors <- c(
  "E.coli" = "#ED7D31",
  "S.cerevisiae" = "#4472C4"
)


ggplot(df, aes(x = species, y = correlation, color = species)) +
  geom_beeswarm(cex = 2, size = 2.5, priority = "density") +
  scale_color_manual(values = species_colors) +
  scale_x_discrete(
    labels = c(
      "E.coli" = expression(italic("E.coli")),
      "S.cerevisiae" = expression(italic("S.cerevisiae"))
    )
  ) +
  labs(
    x = "",
    y = expression("Spearman's " ~ rho),
    title = expression(
      atop(
        "Correlation between " * a[i] * "[" * E[i] * "]"^2 / (a[i+1] * "[" * E[i+1] * "]"^2) * " and " * K[i] *"(EOI)",
        "for glycolysis enzymes in " * italic("E.coli") * " and " * italic("S.cerevisiae")
      )
    )
  ) +
  ylim(0.5, 1) +
  theme(
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20), 
    axis.title.y = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 16),
    plot.title.position = "plot",
    plot.title = element_text(size = 22, hjust = 0.5, vjust = 1),
    panel.grid = element_blank(),  
    panel.background = element_rect(fill = "white", color = NA),  
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    legend.position = "none"   
  )
