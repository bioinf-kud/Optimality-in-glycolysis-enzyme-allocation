library(ggplot2)
library(ggrepel)
library(dplyr)
colour<-c("#4472C4","#ED7D31","#70AD47","#824098","#C77EB5","#FFD166","#B89C6F","#8BBDE3")
dir<-"-path to folder-/figures_and_source_data"
head<-"/figure 4/data/calculated/"
tail1<-".glycolysis.tsv"
selected="Tumor"
data5 <- read.table(paste0(dir,head,"TCGA-COAD",tail1),sep ='\t')
data9 <- read.table(paste0(dir,head,"TCGA-BRCA",tail1),sep='\t')
data13 <- read.table(paste0(dir,head,"TCGA-HNSC",tail1),sep='\t')
data17 <- read.table(paste0(dir,head,"TCGA-KIRC",tail1),sep='\t')
data21 <- read.table(paste0(dir,head,"TCGA-LUAD",tail1),sep='\t')
data25 <- read.table(paste0(dir,head,"TCGA-LUSC",tail1),sep='\t')
data29 <- read.table(paste0(dir,head,"TCGA-PRAD",tail1),sep='\t')
data33 <- read.table(paste0(dir,head,"TCGA-THCA",tail1),sep='\t')
data9<-subset(data9,select=-V10)
data13<-subset(data13,select=-V10)
data17<-subset(data17,select=-V10)
data21<-subset(data21,select=-V10)
data25<-subset(data25,select=-V10)
data29<-subset(data29,select=-V10)
data33<-subset(data33,select=-V10)


data5$cat <- data5$V1
data5<-data5[data5$cat==selected,]
data5<-t(data5)
colnames(data5)<-make.unique(as.character(data5[1, ]))
data5 <- data5[-10, ]
data5 <- data5[-1, ]
data5 <- as.data.frame(data5)
data5[]<-lapply(data5, function(x) as.numeric(as.character(x)))
data9 <- data9[, -10]
data9$V10 <- gsub("\\..*", "",data9$V1)
data9<-data9[data9$V10==selected,]
data9<-subset(data9,select=-V10)
data9 <- t(data9)
colnames(data9) <- data9[1, ]
data9 <- data9[-1, ]
data9 <- as.data.frame(data9)
data9[]<-lapply(data9, function(x) as.numeric(as.character(x)))
data13 <- data13[, -10]
data13$V10 <- gsub("\\..*", "",data13$V1)
data13<-data13[data13$V10==selected,]
data13<-subset(data13,select=-V10)
data13 <- t(data13)
colnames(data13) <- data13[1, ]
data13 <- data13[-1, ]
data13 <- as.data.frame(data13)
data13[]<-lapply(data13, function(x) as.numeric(as.character(x)))
data17 <- data17[, -10]
data17$V10 <- gsub("\\..*", "",data17$V1)
data17<-data17[data17$V10==selected,]
data17<-subset(data17,select=-V10)
data17 <- t(data17)
colnames(data17) <- data17[1, ]
data17 <- data17[-1, ]
data17 <- as.data.frame(data17)
data17[]<-lapply(data17, function(x) as.numeric(as.character(x)))
data21 <- data21[, -10]
data21$V10 <- gsub("\\..*", "",data21$V1)
data21<-data21[data21$V10==selected,]
data21<-subset(data21,select=-V10)
data21 <- t(data21)
colnames(data21) <- data21[1, ]
data21 <- data21[-1, ]
data21 <- as.data.frame(data21)
data21[]<-lapply(data21, function(x) as.numeric(as.character(x)))
data25 <- data25[, -10]
data25$V10 <- gsub("\\..*", "",data25$V1)
data25<-data25[data25$V10==selected,]
data25<-subset(data25,select=-V10)
data25 <- t(data25)
colnames(data25) <- data25[1, ]
data25 <- data25[-1, ]
data25 <- as.data.frame(data25)
data25[]<-lapply(data25, function(x) as.numeric(as.character(x)))
data29 <- data29[, -10]
data29$V10 <- gsub("\\..*", "",data29$V1)
data29<-data29[data29$V10==selected,]
data29<-subset(data29,select=-V10)
data29 <- t(data29)
colnames(data29) <- data29[1, ]
data29 <- data29[-1, ]
data29 <- as.data.frame(data29)
data29[]<-lapply(data29, function(x) as.numeric(as.character(x)))
data33 <- data33[, -10]
data33$V10 <- gsub("\\..*", "",data33$V1)
data33<-data33[data33$V10==selected,]
data33<-subset(data33,select=-V10)
data33 <- t(data33)
colnames(data33) <- data33[1, ]
data33 <- data33[-1, ]
data33 <- as.data.frame(data33)
data33[]<-lapply(data33, function(x) as.numeric(as.character(x)))
vector <- c(8.841552744, -0.6857715, 7.833065244, -9.641140506, -6.172457509, 5.478969765, -1.774938001, 1.290864001)





data5 <- log(data5)
data <- data5 %>%
  rowwise() %>%
  summarise(
    mean = mean(c_across(everything())),  
    sd = sd(c_across(everything()))     
  )
data<-as.data.frame(data)
data$K <- vector
data$pair <- c("Glc->G6P->F6P", "G6P->F6P->FBP", "F6P->FBP->DHAP+GAP", "FBP->GAP->BPG", "GAP->BPG->3PG", "BPG->3PG->2PG", "3PG->2PG->PEP", "2PG->PEP->Pyr")
data$samp="COAD"
summary<-data
data9 <- log(data9)
data <- data9 %>%
  rowwise() %>%
  summarise(
    mean = mean(c_across(everything())),  
    sd = sd(c_across(everything()))     
  )
data<-as.data.frame(data)
data$K <- vector
data$pair <- c("Glc->G6P->F6P", "G6P->F6P->FBP", "F6P->FBP->DHAP+GAP", "FBP->GAP->BPG", "GAP->BPG->3PG", "BPG->3PG->2PG", "3PG->2PG->PEP", "2PG->PEP->Pyr")
data$samp="BRCA"
summary<-rbind(summary,data)
data13 <- log(data13)
data <- data13 %>%
  rowwise() %>%
  summarise(
    mean = mean(c_across(everything())),  
    sd = sd(c_across(everything()))     
  )
data<-as.data.frame(data)
data$K <- vector
data$pair <- c("Glc->G6P->F6P", "G6P->F6P->FBP", "F6P->FBP->DHAP+GAP", "FBP->GAP->BPG", "GAP->BPG->3PG", "BPG->3PG->2PG", "3PG->2PG->PEP", "2PG->PEP->Pyr")
data$samp="HNSC"
summary<-rbind(summary,data)
data17 <- log(data17)
data <- data17 %>%
  rowwise() %>%
  summarise(
    mean = mean(c_across(everything())),  
    sd = sd(c_across(everything()))     
  )
data<-as.data.frame(data)
data$K <- vector
data$pair <- c("Glc->G6P->F6P", "G6P->F6P->FBP", "F6P->FBP->DHAP+GAP", "FBP->GAP->BPG", "GAP->BPG->3PG", "BPG->3PG->2PG", "3PG->2PG->PEP", "2PG->PEP->Pyr")
data$samp="KIRC"
summary<-rbind(summary,data)
data21 <- log(data21)
data <- data21 %>%
  rowwise() %>%
  summarise(
    mean = mean(c_across(everything())),  
    sd = sd(c_across(everything()))     
  )
data<-as.data.frame(data)
data$K <- vector
data$pair <- c("Glc->G6P->F6P", "G6P->F6P->FBP", "F6P->FBP->DHAP+GAP", "FBP->GAP->BPG", "GAP->BPG->3PG", "BPG->3PG->2PG", "3PG->2PG->PEP", "2PG->PEP->Pyr")
data$samp="LUAD"
summary<-rbind(summary,data)
data25 <- log(data25)
data <- data25 %>%
  rowwise() %>%
  summarise(
    mean = mean(c_across(everything())),  
    sd = sd(c_across(everything()))     
  )
data<-as.data.frame(data)
data$K <- vector
data$pair <- c("Glc->G6P->F6P", "G6P->F6P->FBP", "F6P->FBP->DHAP+GAP", "FBP->GAP->BPG", "GAP->BPG->3PG", "BPG->3PG->2PG", "3PG->2PG->PEP", "2PG->PEP->Pyr")
data$samp="LUSC"
summary<-rbind(summary,data)
data29 <- log(data29)
data <- data29 %>%
  rowwise() %>%
  summarise(
    mean = mean(c_across(everything())),  
    sd = sd(c_across(everything()))     
  )
data<-as.data.frame(data)
data$K <- vector
data$pair <- c("Glc->G6P->F6P", "G6P->F6P->FBP", "F6P->FBP->DHAP+GAP", "FBP->GAP->BPG", "GAP->BPG->3PG", "BPG->3PG->2PG", "3PG->2PG->PEP", "2PG->PEP->Pyr")
data$samp="PRAD"
summary<-rbind(summary,data)
data33 <- log(data33)
data <- data33 %>%
  rowwise() %>%
  summarise(
    mean = mean(c_across(everything())),  
    sd = sd(c_across(everything()))     
  )
data<-as.data.frame(data)
data$K <- vector
data$pair <- c("Glc->G6P->F6P", "G6P->F6P->FBP", "F6P->FBP->DHAP+GAP", "FBP->GAP->BPG", "GAP->BPG->3PG", "BPG->3PG->2PG", "3PG->2PG->PEP", "2PG->PEP->Pyr")
data$samp="THCA"
summary<-rbind(summary,data)

ggplot(summary) + 
  geom_point(aes(y = K, x = mean, color = samp), size = 8, alpha = 0.5) +  
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#ED7D31", alpha = 0.5) + 
  labs(
    title = expression(
      atop(
        "Comparison of mean of " * a[i] * "[" * E[i] * "]"^2 / (a[i+1] * "[" * E[i+1] * "]"^2) * " and " * K[i],
        "for glycolysis enzymes in human tumor samples"
      )
    ),
    x = expression("Log"* (a[i] * "[" * E[i] * "]"^2 / (a[i+1] * "[" * E[i+1] * "]"^2))),  
    y = expression("Log" *(K[i])),
    col = ""
  ) +
  scale_color_manual(
    values = colour,
    labels = c("BRCA (n=1117)", "COAD (n=291)", "HNSC (n=521)", "KIRC (n=541)", 
               "LUAD (n=540)", "LUSC (n=501)", "PRAD (n=501)", "THCA (n=512)")  
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
    legend.position = "right", 
    legend.background = element_rect(fill = "white", color = NA), 
    legend.text = element_text(size = 16)  
  ) +
  guides(fill = guide_legend(title = NULL))  +
  xlim(-15, 15) +
  ylim(-15, 15) +
  geom_text_repel(
    data = subset(summary, summary$samp == "THCA"),
    aes(label = pair), size = 7, 
    y = summary[summary$samp == "THCA", ]$K, 
    x = summary[summary$samp == "THCA", ]$mean
  )

