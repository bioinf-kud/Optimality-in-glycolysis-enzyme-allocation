library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
vector <- c(8.841552744, -0.6857715, 7.833065244, -9.641140506, -6.172457509, 5.478969765, -1.774938001, 1.290864001)
dir<-"-path to folder-/figures_and_source_data"
head<-"/figure 4/data/calculated/CCLE_proteome.glycolysis"
tail<-".tsv"
data1 <- read.table(paste0(dir,head,tail),header=FALSE,sep='\t')
row.names(data1)<-data1$V1
data1<-subset(data1,select=-V1)
data1<-subset(data1,select=-V10)
data1<-t(data1)
data1<-as.data.frame(data1)
data1 <- log(data1)
data <- data1 %>%
  rowwise() %>%
  summarise(
    mean = mean(c_across(everything())),  
    sd = sd(c_across(everything()))     
  )
data<-as.data.frame(data)
data$K <- vector
data$pair <- c("Glc->G6P->F6P", "G6P->F6P->FBP", "F6P->FBP->DHAP+GAP", "FBP->GAP->BPG", "GAP->BPG->3PG", "BPG->3PG->2PG", "3PG->2PG->PEP", "2PG->PEP->Pyr")


ggplot(data) + 
  geom_point(aes(y = K, x = mean), size = 10,alpha = 0.5, color = "#ED7D31") +  
  geom_errorbarh(aes(xmin = mean - sd, xmax = mean + sd, y = K), height = 0.5, color = "#ED7D31")+
  geom_text_repel(
    data = data,
    aes(label = pair), size = 7,y =data$K, x = data$mean)+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#ED7D31", alpha = 0.5) + 
  labs(
    title = expression(
      atop(
        "Comparison of " * a[i] * "[" * E[i] * "]"^2 / (a[i+1] * "[" * E[i+1] * "]"^2) * " and " * K[i],
        "for glycolysis enzymes in 378 human cell lines"
      )
    ),
    x = expression("Log"* (a[i] * "[" * E[i] * "]"^2 / (a[i+1] * "[" * E[i+1] * "]"^2))),  
    y = expression("Log" *(K[i])),
    col = "method"
  ) +
  theme(
    axis.text.x = element_text(size = 20),  
    axis.title.x = element_text(size = 20),  
    axis.text.y = element_text(size = 20),  
    axis.title = element_text(size = 20),  
    plot.title = element_text(size = 22, hjust = 0.5, vjust = 1), 
    plot.title.position = "plot",  
    panel.grid = element_blank(),  
    panel.background = element_rect(fill = "white", color = NA),  
    panel.border = element_rect(color = "black", fill = NA, size = 1),  
    legend.position = c(0.85, 0.17), 
    legend.background = element_rect(fill = "white", color = NA), 
    legend.text = element_text(size = 16)  
  ) +
  guides(fill = guide_legend(title = NULL), colour = "none")  +
  xlim(-15, 15) +
  ylim(-15, 15)
           
