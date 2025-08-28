library(ggplot2)
library(ggrepel)
dir<-"-path to folder-/figures_and_source_data"
head<-"/figureS6/data/calculated/eco.glycolysis"
tail<-".tsv"
data1 <- read.table(paste0(dir,head,tail),sep='\t')
data2 <- read.table(paste0(dir,head,".a_permuted_100",tail),sep='\t')
data3 <- read.table(paste0(dir,head,".e_permuted_100",tail),sep='\t')
data4 <- read.table(paste0(dir,head,".both_permuted_100",tail),sep='\t')
row.names(data1)<-data1$V1
data1<-subset(data1,select=-V1)
data1<-t(data1)
data1<-as.data.frame(data1)

data2<-subset(data2,select=-V1)
data2<-t(data2)
data2<-as.data.frame(data2)
data3<-subset(data3,select=-V1)
data3<-t(data3)
data3<-as.data.frame(data3)
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
correlations$datanum <-datanum
correlations$per<- "original data"
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
correlations$datanum <-datanum
correlations$per<-"a permuted"
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
correlations$datanum <-datanum
correlations$per<-"[E] permuted"
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
correlations$datanum <-datanum
correlations$per<-"both permuted"
Bmean=mean(na.omit(correlations$value))
dataf<-rbind(dataf,correlations)



ggplot(dataf) + 
  geom_boxplot(aes(y = value, x = per, fill = per, colour = per), alpha = 0.5) +  
  geom_jitter(aes(y = value, x = per, colour = per), size = 1, alpha = 0.1) +  
  scale_fill_manual(
    values = c(
      "both permuted" = "#4472C4", 
      "[E] permuted" = "#ED7D31", 
      "a permuted" = "#70AD47", 
      "original data" = "#824098"
    )
  ) +
  scale_colour_manual(
    values = c(
      "both permuted" = "#4472C4", 
      "[E] permuted" = "#ED7D31", 
      "a permuted" = "#70AD47", 
      "original data" = "#824098"
    )
  ) +
  labs(
    title = expression(
      atop(
        "Correlation between " * a[i] * "[" * E[i] * "]"^2 / (a[i+1] * "[" * E[i+1] * "]"^2) * " and " * K[i],
        "for glycolysis enzymes in " * italic("E.coli")
      )
    ),
    y = expression("Spearman's " ~ rho),  
    col = "method"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),  
    axis.title.x = element_blank(),  
    axis.text.y = element_text(size = 20),  
    axis.title = element_text(size = 20),  
    plot.title = element_text(size = 22, hjust = 0.5, vjust = 1), 
    plot.title.position = "plot",  
    panel.grid = element_blank(),  
    panel.background = element_rect(fill = "white", color = NA),  
    panel.border = element_rect(color = "black", fill = NA, size = 1),  
    legend.position = "none",
  ) +
  guides(fill = guide_legend(title = NULL), colour = "none")  
