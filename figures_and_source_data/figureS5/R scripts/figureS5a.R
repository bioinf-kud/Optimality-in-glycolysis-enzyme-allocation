library(ggplot2)
library(ggrepel)
dir<-"-path to folder-/figures_and_source_data"
head<-"/figureS5/data/calculated/"
data <- read.table(paste0(dir,head,"CCLE_proteome.glycolysis.tsv"), sep = "\t")
data1 <- read.table(paste0(dir,head,"CCLE_proteome.glycolysis.no_substrate.tsv"), sep = "\t")
vector <- c(8.841552744, -0.6857715, 7.833065244, -9.641140506, -6.172457509, 5.478969765, -1.774938001, 1.290864001)
colour<-c("#4472C4","#ED7D31")
row.names(data)<-data$V1
row.names(data1)<-data1$V1
data<-subset(data,select=-V1)
data<-subset(data,select=-V10)
data1<-subset(data1,select=-V1)
data1<-subset(data1,select=-V10)
data<-t(data)
data1<-t(data1)
data<-as.data.frame(data)
data1<-as.data.frame(data1)
log_data <- log(data)
log_data1 <- log(data1) 
correlations <- numeric(ncol(log_data))
correlations1 <- numeric(ncol(log_data1))
for (i in 1:ncol(log_data)) {
  correlations[i] <- cor(log_data[[i]], vector, method = "spearman")
}
correlations<-as.data.frame(correlations)
row.names(correlations)<-names(data)
for (i in 1:ncol(log_data1)) {
  correlations1[i] <- cor(log_data1[[i]], vector, method = "spearman")
}
correlations1<-as.data.frame(correlations1)
row.names(correlations1)<-names(data1)
correlations$method="with_substrate"
Amean<-mean(correlations$correlations)
correlations1$method="without_substrate"
Bmean<-mean(correlations1$correlations)
names(correlations1)<-c("correlations","method")
correlations<-rbind(correlations,correlations1)
ggplot(correlations) + 
  geom_density(aes(x = correlations, fill = method), colour = NA) + 
  ggtitle("Spearman's correlations distribution(n=378)") +
  xlim(-1,1)

ggplot(correlations) +
  geom_histogram(aes(x = correlations, fill = method, colour = method), binwidth = 0.025, position = "identity", alpha = 0.5) +  # 不堆积，透明度为 0.5
  scale_fill_manual(
    values = c("#4472C4","#ED7D31"),  # 填充颜色
    labels = c(
      bquote("Adjusted a (Spearman's " ~ bar(rho) ~ ":" ~ .(sprintf("%.2f", Amean)) ~ ")"),
      bquote("Original a (Spearman's " ~ bar(rho) ~ ":" ~ .(sprintf("%.2f", Bmean)) ~ ")")
    )
  )+
  scale_colour_manual(values = c("#4472C4","#ED7D31")) +
  labs(
    title = expression(
      atop(
        "Correlation between " * a[i] * "[" * E[i] * "]"^2 / (a[i+1] * "[" * E[i+1] * "]"^2) * " and " * K[i],
        "for glycolysis enzymes in 378 human cell lines"
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

