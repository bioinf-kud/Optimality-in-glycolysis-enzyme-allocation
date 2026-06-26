library(ggplot2)
library(ggrepel)
library(ggbeeswarm)
library(dplyr)
library(tidyr)
get_pearson <- function(data, vector,spec) {
  correlations <- numeric(ncol(data))
  for (i in 1:ncol(data)) {
    current_column <- as.numeric(data[[i]])
    current_column <- log(current_column)
    valid_indices <- which(current_column != -Inf & current_column != Inf)
    filtered_column <- current_column[valid_indices]
    filtered_vector <- vector[valid_indices]
    if (length(valid_indices) > 1) {
      correlations[i] <- cor(filtered_column, filtered_vector, method = "pearson")
    } else {
      correlations[i] <- NA
    }
  }
  data.frame(
    correlation_Pearson = correlations,
    species = spec
  )
}

vector <- c(13.34312772,
            -0.077827744,
            8.924192273,
            -0.29214456,
            -2.006169471,
            -3.01912952,
            3.362797029,
            -1.682640295,
            1.171297752
)
vector <- vector / log(10)
dir<-"-path to folder-/figures_and_source_data"
head<-"/figure 4/data/E.coli_calc"
tail<-".csv"
data1 <- read.csv(paste0(dir,head,tail),header=TRUE)
eco<-get_pearson(data1, vector,"E.coli")

vector <- c(5.522477554, 8.111677106,-0.6857715,7.103189605,-0.29214456,-3.065802002,-5.00598449,6.209460215,-1.774938001,1.290864001)
vector <- vector / log(10)
head<-"/figure 4/data/S.cerevisiae_calc"
data1 <- read.csv(paste0(dir,head,tail),header=TRUE)
sce<-get_pearson(data1, vector,"S.cerevisiae")
df<-rbind(eco,sce)


species_colors <- c(
  "E.coli" = "#ED7D31",
  "S.cerevisiae" = "#4472C4"
)


p<-ggplot(df, aes(x = species, y = correlation_Pearson, color = species)) +
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
    y = expression("Pearson's r"),
    title = expression(
      atop(
        "Correlation between K and CAQ (EOI) for",
        "glycolysis enzymes in " * italic("E.coli") * " and " * italic("S.cerevisiae")
      )
    )
  ) +
  ylim(-0.2, 1) +
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
ggsave(paste0(dir,"/figureS5/figures/","figureS5b.pdf"), p, width = 650/96, height = 700/96, units = "in", device = "pdf")
