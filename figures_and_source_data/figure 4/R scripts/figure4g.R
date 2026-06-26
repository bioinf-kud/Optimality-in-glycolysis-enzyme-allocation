library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
vector <- c(5.522477554, 8.111677106,-0.6857715,7.103189605,-0.29214456,-3.065802002,-5.00598449,6.209460215,-1.774938001,1.290864001)
vector <- vector / log(10)
dir<-"-path to folder-/figures_and_source_data"
head<-"/figure 4/data/S.cerevisiae_calc"
tail<-".csv"
data1 <- read.csv(paste0(dir,head,tail),header=TRUE)
data1<-log(data1,10)
data <- data1 %>%
  rowwise() %>%
  summarise(
    mean = mean(c_across(everything())),  
    sd = sd(c_across(everything()))     
  )
data<-as.data.frame(data)
data$K <- vector
data$pair <- c("GLUT-HK","HK-GPI", "GPI-PFK", "PFK-ALDO", "ALDO-TPI", "TPI-GAPDH", "GAPDH-PGK", "PGK-PGAM", "PGAM-ENO", "ENO-PK")

line_df <- data.frame(
  intercept = c(0, 1, -1, 2, -2),
  level = factor(
    c(
      "Equal",
      "10-fold difference",
      "10-fold difference",
      "100-fold difference",
      "100-fold difference"
    ),
    levels = c(
      "Equal",
      "10-fold difference",
      "100-fold difference"
    )
  )
)

p<-ggplot(data) + 
  geom_point(aes(y = K, x = mean), size = 10, alpha = 0.5, color = "#ED7D31") +  
  geom_errorbarh(aes(xmin = mean - sd, xmax = mean + sd, y = K), height = 0.5, color = "#ED7D31") +
  geom_text_repel(
    data = data,
    aes(label = pair), size = 7, y = data$K, x = data$mean
  ) +
  geom_abline(
    data = line_df,
    aes(slope = 1, intercept = intercept, alpha = level),
    linetype = "dashed",
    color = "#4472C4"
  ) +
  scale_alpha_manual(
    name = "",
    values = c(
      "Equal" = 1,
      "10-fold difference" = 0.7,
      "100-fold difference" = 0.3
    )
  ) +
  guides(alpha = guide_legend(
    override.aes = list(linetype = "dashed", color = "#4472C4")
  )) +
  labs(
    title = expression(
      atop(
        "Comparison of K and CAQ for glycolysis" ,
        "enzymes in 16 "* italic("S.cerevisiae")*" samples"
      )
    ),
    x = expression(log[10]* (CAQ)),  
    y = expression(log[10] *(K)),
    col = "method"
  ) +
  theme(
    axis.text.x = element_text(size = 20,colour = "black"),  
    axis.title.x = element_text(size = 20,colour = "black"),  
    axis.text.y = element_text(size = 20,colour = "black"),  
    axis.title = element_text(size = 20,colour = "black"),  
    axis.ticks = element_line(colour = "black"),
    axis.line = element_line(colour = "black", linewidth = 0.8),
    plot.title = element_text(size = 20, hjust = 0.5, vjust = 1), 
    plot.title.position = "plot",  
    panel.grid = element_blank(),  
    panel.background = element_rect(fill = "white", color = NA),  
    panel.border = element_rect(color = "black", fill = NA, size = 1),  
    legend.position = c(0.72, 0.11), 
    legend.background = element_rect(fill = "white", color = NA), 
    legend.text = element_text(size = 20)
  ) +
  xlim(-6, 6) +
  ylim(-6, 6)
ggsave(paste0(dir,"/figure 4/figures/","figure 4g.pdf"), p, width = 650/96, height = 700/96, units = "in", device = "pdf")


