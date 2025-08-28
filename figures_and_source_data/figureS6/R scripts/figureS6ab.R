library(ggplot2)
library(ggrepel)
library(ggbeeswarm)
dir<-"-path to folder-/figures_and_source_data"
head<-"/figureS6/data/"
tail<-".csv"
df <- read.csv(paste0(dir,head,"corr",tail),header=TRUE)
df1<- read.csv(paste0(dir,head,"correlation_without_substrate",tail),header=TRUE)
df$method="with_substrate"
Amean<-mean(df[df$species=="E.coli",]$correlation)
Bmean<-mean(df[df$species=="S.cerevisiae",]$correlation)
df1$method="without_substrate"
Cmean<-mean(df1[df1$species=="E.coli",]$correlation)
Dmean<-mean(df1[df1$species=="S.cerevisiae",]$correlation)
df<-rbind(df,df1)
ggplot(correlations) + 
  geom_density(aes(x = correlations, fill = method), colour = NA) + 
  ggtitle("Spearman's correlations distribution(n=378)") +
  xlim(-1,1)
correlations<-df[df$species=="E.coli",]
ggplot(correlations) +
  geom_histogram(aes(x = correlation, fill = method, colour = method), binwidth = 0.025, position = "identity", alpha = 0.5) +  # 不堆积，透明度为 0.5
  scale_fill_manual(
    values = c("#4472C4","#ED7D31"),  # 填充颜色
    labels = c(
      bquote("Adjusted a (Spearman's " ~ bar(rho) ~ ":" ~ .(sprintf("%.2f", Amean)) ~ ")"),
      bquote("Original a (Spearman's " ~ bar(rho) ~ ":" ~ .(sprintf("%.2f", Cmean)) ~ ")")
    )
  )+
  scale_colour_manual(values = c("#4472C4","#ED7D31")) +
  labs(
    title = expression(
      atop(
        "Correlation between " * a[i] * "[" * E[i] * "]"^2 / (a[i+1] * "[" * E[i+1] * "]"^2) * " and " * K[i],
        "for glycolysis enzymes in 66 "* italic("E.coli")*" samples"
      )
    ),
    x = expression("Spearman's " ~ rho),
    y = "Number of cell lines",
    col = "method"
  ) +
  xlim(-1, 1) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 22, hjust = 0.5, vjust = 1),
    plot.title.position = "plot",
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    legend.position = c(0.37, 0.87),  # 图例位置调整到左上角
    legend.background = element_rect(fill = "white", color = NA),  # 图例背景设置为白色
    legend.text = element_text(size = 16)  # 图例字体大小
  ) +
  guides(fill = guide_legend(title = NULL), colour = "none")  # 移除图例标题






correlations<-df[df$species=="S.cerevisiae",]
ggplot(correlations) +
  geom_histogram(aes(x = correlation, fill = method, colour = method), binwidth = 0.025, position = "identity", alpha = 0.5) +  # 不堆积，透明度为 0.5
  scale_fill_manual(
    values = c("#4472C4","#ED7D31"),  # 填充颜色
    labels = c(
      bquote("Adjusted a (Spearman's " ~ bar(rho) ~ ":" ~ .(sprintf("%.2f", Bmean)) ~ ")"),
      bquote("Original a (Spearman's " ~ bar(rho) ~ ":" ~ .(sprintf("%.2f", Dmean)) ~ ")")
    )
  )+
  scale_colour_manual(values = c("#4472C4","#ED7D31")) +
  labs(
    title = expression(
      atop(
        "Correlation between " * a[i] * "[" * E[i] * "]"^2 / (a[i+1] * "[" * E[i+1] * "]"^2) * " and " * K[i],
        "for glycolysis enzymes in 16 "* italic("S.cerevisiae")*" samples"
      )
    ),
    x = expression("Spearman's " ~ rho),
    y = "Number of cell lines",
    col = "method"
  ) +
  xlim(-1, 1) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 22, hjust = 0.5, vjust = 1),
    plot.title.position = "plot",
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    legend.position = c(0.37, 0.87),  # 图例位置调整到左上角
    legend.background = element_rect(fill = "white", color = NA),  # 图例背景设置为白色
    legend.text = element_text(size = 16)  # 图例字体大小
  ) +
  guides(fill = guide_legend(title = NULL), colour = "none")  # 移除图例标题

