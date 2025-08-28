library(ggplot2)
library(ggrepel)
dir<-"-path to folder-/figures_and_source_data"
head<-"/figureS5/data/calculated/CCLE_proteome.glycolysis"
tail<-".tsv"
data1 <- read.table(paste0(dir,head,tail),sep='\t')
data2 <- read.table(paste0(dir,head,".a_permuted",tail),sep='\t')
data3 <- read.table(paste0(dir,head,".e_permuted",tail),sep='\t')
data4 <- read.table(paste0(dir,head,".both_permuted",tail),sep='\t')
row.names(data1)<-data1$V1
data1<-subset(data1,select=-V1)
data1<-t(data1)
data1<-as.data.frame(data1)
row.names(data2)<-data1$V1
data2<-subset(data2,select=-V1)
data2<-t(data2)
data2<-as.data.frame(data2)
row.names(data3)<-data3$V1
data3<-subset(data3,select=-V1)
data3<-t(data3)
data3<-as.data.frame(data3)
row.names(data4)<-data4$V1
data4<-subset(data4,select=-V1)
data4<-t(data4)
data4<-as.data.frame(data4)
vector <- c(8.841552744, -0.6857715, 7.833065244, -9.641140506, -6.172457509, 5.478969765, -1.774938001, 1.290864001)
data<-data1
correlations <- numeric(ncol(data))
datanum<- numeric(ncol(data))
for (i in 1:ncol(data)) {
  current_column <- data[[i]]
  current_column<- log(current_column)
  valid_indices <- which(current_column != -Inf)
  filtered_column <- current_column[valid_indices]
  filtered_vector <- vector[valid_indices]
  valid_indices <- which(current_column != Inf)
  filtered_column <- filtered_column[valid_indices]
  filtered_vector <- filtered_vector[valid_indices]
  if (length(valid_indices) > 1) {
    datanum[i]<- length(valid_indices)
    correlations[i] <- cor(filtered_column, filtered_vector, method = "spearman")
  } else {
    correlations[i] <- NA
  }
}
correlations <- as.data.frame(correlations)
names(correlations)<-c("value")
row.names(correlations) <- names(data)
correlations$datanum <-datanum
correlations$per<- "Original data"
Omean=mean(correlations$value)
dataf<-correlations
data<-data2
correlations <- numeric(ncol(data))
datanum<- numeric(ncol(data))
for (i in 1:ncol(data)) {
  current_column <- data[[i]]
  current_column<- log(current_column)
  valid_indices <- which(current_column != -Inf)
  filtered_column <- current_column[valid_indices]
  filtered_vector <- vector[valid_indices]
  valid_indices <- which(current_column != Inf)
  filtered_column <- filtered_column[valid_indices]
  filtered_vector <- filtered_vector[valid_indices]
  if (length(valid_indices) > 1) {
    datanum[i]<- length(valid_indices)
    correlations[i] <- cor(filtered_column, filtered_vector, method = "spearman")
  } else {
    correlations[i] <- NA
  }
}
correlations <- as.data.frame(correlations)
names(correlations)<-c("value")
row.names(correlations) <- names(data)
correlations$datanum <-datanum
correlations$per<-"Enzyme kinetics permuted"
Amean=mean(correlations$value)
dataf<-rbind(dataf,correlations)
data<-data3
correlations <- numeric(ncol(data))
datanum<- numeric(ncol(data))
for (i in 1:ncol(data)) {
  current_column <- data[[i]]
  current_column<- log(current_column)
  valid_indices <- which(current_column != -Inf)
  filtered_column <- current_column[valid_indices]
  filtered_vector <- vector[valid_indices]
  valid_indices <- which(current_column != Inf)
  filtered_column <- filtered_column[valid_indices]
  filtered_vector <- filtered_vector[valid_indices]
  if (length(valid_indices) > 1) {
    datanum[i]<- length(valid_indices)
    correlations[i] <- cor(filtered_column, filtered_vector, method = "spearman")
  } else {
    correlations[i] <- NA
  }
}
correlations <- as.data.frame(correlations)
names(correlations)<-c("value")
row.names(correlations) <- names(data)
correlations$datanum <-datanum
correlations$per<-"Enzyme concentration permuted"
Emean=mean(na.omit(correlations$value))
dataf<-rbind(dataf,correlations)
data<-data4
correlations <- numeric(ncol(data))
datanum<- numeric(ncol(data))
for (i in 1:ncol(data)) {
  current_column <- data[[i]]
  current_column<- log(current_column)
  valid_indices <- which(current_column != -Inf)
  filtered_column <- current_column[valid_indices]
  filtered_vector <- vector[valid_indices]
  valid_indices <- which(current_column != Inf)
  filtered_column <- filtered_column[valid_indices]
  filtered_vector <- filtered_vector[valid_indices]
  if (length(valid_indices) > 1) {
    datanum[i]<- length(valid_indices)
    correlations[i] <- cor(filtered_column, filtered_vector, method = "spearman")
  } else {
    correlations[i] <- NA
  }
}
correlations <- as.data.frame(correlations)
names(correlations)<-c("value")
row.names(correlations) <- names(data)
correlations$datanum <-datanum
correlations$per<-"Both permuted"
Bmean=mean(na.omit(correlations$value))
dataf<-rbind(dataf,correlations)
dataf

ggplot(dataf) +
  geom_histogram(aes(x = value, fill = per, colour = per), binwidth = 0.025, position = "identity", alpha = 0.5) +  # 不堆积，透明度为 0.5
  scale_fill_manual(
    values = c("#4472C4","#ED7D31","#70AD47","#824098"),  
    labels = c(
      bquote("both permuted (Spearman's " ~ bar(rho) ~ ":" ~ .(sprintf("%.2f", Bmean)) ~ ")"),
      bquote("[E] permuted (Spearman's " ~ bar(rho) ~ ":" ~ .(sprintf("%.2f", Emean)) ~ ")"),
      bquote("a permuted (Spearman's " ~ bar(rho) ~ ":" ~ .(sprintf("%.2f", Amean)) ~ ")"),
      bquote("original data (Spearman's " ~ bar(rho) ~ ":" ~ .(sprintf("%.2f", Omean)) ~ ")")
    )
  )+
  scale_colour_manual(values = c("#4472C4","#ED7D31","#70AD47","#824098")) +
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
    legend.position = c(0.36, 0.87),  
    legend.background = element_rect(fill = "transparent", color = NA), 
    legend.text = element_text(size = 16)  
  ) +
  guides(fill = guide_legend(title = NULL), colour = "none")  
