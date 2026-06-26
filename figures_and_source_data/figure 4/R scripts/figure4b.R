library(ggplot2)
library(ggrepel)
dir<-"-path to folder-/figures_and_source_data"
head<-"/figure 4/data/calculated/"
tail1<-".glycolysis.tsv"
selected="Normal"
colour<-c("#4472C4","#ED7D31","#70AD47","#824098","#C77EB5","#FFD166","#B89C6F","#8BBDE3")
data5 <- read.table(paste0(dir,head,"TCGA-COAD",tail1),sep ='\t')
data9 <- read.table(paste0(dir,head,"TCGA-BRCA",tail1),sep='\t')
data13 <- read.table(paste0(dir,head,"TCGA-HNSC",tail1),sep='\t')
data17 <- read.table(paste0(dir,head,"TCGA-KIRC",tail1),sep='\t')
data21 <- read.table(paste0(dir,head,"TCGA-LUAD",tail1),sep='\t')
data25 <- read.table(paste0(dir,head,"TCGA-LUSC",tail1),sep='\t')
data29 <- read.table(paste0(dir,head,"TCGA-PRAD",tail1),sep='\t')
data33 <- read.table(paste0(dir,head,"TCGA-THCA",tail1),sep='\t')
data9<-subset(data9,select=-V12)
data13<-subset(data13,select=-V12)
data17<-subset(data17,select=-V12)
data21<-subset(data21,select=-V12)
data25<-subset(data25,select=-V12)
data29<-subset(data29,select=-V12)
data33<-subset(data33,select=-V12)


data5$cat <- data5$V1
data5<-data5[data5$cat==selected,]
data5<-t(data5)
colnames(data5)<-make.unique(as.character(data5[1, ]))
data5 <- data5[-12, ]
data5 <- data5[-1, ]
data5 <- as.data.frame(data5)
data5[]<-lapply(data5, function(x) as.numeric(as.character(x)))
data9 <- data9[, -12]
data9$V12 <- gsub("\\..*", "",data9$V1)
data9<-data9[data9$V12==selected,]
data9<-subset(data9,select=-V12)
data9 <- t(data9)
colnames(data9) <- data9[1, ]
data9 <- data9[-1, ]
data9 <- as.data.frame(data9)
data9[]<-lapply(data9, function(x) as.numeric(as.character(x)))
data13 <- data13[, -12]
data13$V12 <- gsub("\\..*", "",data13$V1)
data13<-data13[data13$V12==selected,]
data13<-subset(data13,select=-V12)
data13 <- t(data13)
colnames(data13) <- data13[1, ]
data13 <- data13[-1, ]
data13 <- as.data.frame(data13)
data13[]<-lapply(data13, function(x) as.numeric(as.character(x)))
data17 <- data17[, -12]
data17$V12 <- gsub("\\..*", "",data17$V1)
data17<-data17[data17$V12==selected,]
data17<-subset(data17,select=-V12)
data17 <- t(data17)
colnames(data17) <- data17[1, ]
data17 <- data17[-1, ]
data17 <- as.data.frame(data17)
data17[]<-lapply(data17, function(x) as.numeric(as.character(x)))
data21 <- data21[, -12]
data21$V12 <- gsub("\\..*", "",data21$V1)
data21<-data21[data21$V12==selected,]
data21<-subset(data21,select=-V12)
data21 <- t(data21)
colnames(data21) <- data21[1, ]
data21 <- data21[-1, ]
data21 <- as.data.frame(data21)
data21[]<-lapply(data21, function(x) as.numeric(as.character(x)))
data25 <- data25[, -12]
data25$V12 <- gsub("\\..*", "",data25$V1)
data25<-data25[data25$V12==selected,]
data25<-subset(data25,select=-V12)
data25 <- t(data25)
colnames(data25) <- data25[1, ]
data25 <- data25[-1, ]
data25 <- as.data.frame(data25)
data25[]<-lapply(data25, function(x) as.numeric(as.character(x)))
data29 <- data29[, -12]
data29$V12 <- gsub("\\..*", "",data29$V1)
data29<-data29[data29$V12==selected,]
data29<-subset(data29,select=-V12)
data29 <- t(data29)
colnames(data29) <- data29[1, ]
data29 <- data29[-1, ]
data29 <- as.data.frame(data29)
data29[]<-lapply(data29, function(x) as.numeric(as.character(x)))
data33 <- data33[, -12]
data33$V12 <- gsub("\\..*", "",data33$V1)
data33<-data33[data33$V12==selected,]
data33<-subset(data33,select=-V12)
data33 <- t(data33)
colnames(data33) <- data33[1, ]
data33 <- data33[-1, ]
data33 <- as.data.frame(data33)
data33[]<-lapply(data33, function(x) as.numeric(as.character(x)))
vector <- c(5.522477554, 8.841552744, -0.6857715, 7.833065244, -0.772381588,
            -3.06580201, -6.172457509, 5.478969765, -1.774938001, 1.290864001)
vector <- vector / log(10)




data5 <- log(data5,10)
data <- data5 %>%
  rowwise() %>%
  summarise(
    mean = mean(c_across(everything())),  
    sd = sd(c_across(everything()))     
  )
data<-as.data.frame(data)
data$K <- vector
data$pair <- c("GLUT-HK","HK-GPI", "GPI-PFK", "PFK-ALDO", "ALDO-TPI", "TPI-GAPDH", "GAPDH-PGK", "PGK-PGAM", "PGAM-ENO", "ENO-PK")
data$samp="COAD"
summary<-data
data9 <- log(data9,10)
data <- data9 %>%
  rowwise() %>%
  summarise(
    mean = mean(c_across(everything())),  
    sd = sd(c_across(everything()))     
  )
data<-as.data.frame(data)
data$K <- vector
data$pair <- c("GLUT-HK","HK-GPI", "GPI-PFK", "PFK-ALDO", "ALDO-TPI", "TPI-GAPDH", "GAPDH-PGK", "PGK-PGAM", "PGAM-ENO", "ENO-PK")
data$samp="BRCA"
summary<-rbind(summary,data)
data13 <- log(data13,10)
data <- data13 %>%
  rowwise() %>%
  summarise(
    mean = mean(c_across(everything())),  
    sd = sd(c_across(everything()))     
  )
data<-as.data.frame(data)
data$K <- vector
data$pair <- c("GLUT-HK","HK-GPI", "GPI-PFK", "PFK-ALDO", "ALDO-TPI", "TPI-GAPDH", "GAPDH-PGK", "PGK-PGAM", "PGAM-ENO", "ENO-PK")
data$samp="HNSC"
summary<-rbind(summary,data)
data17 <- log(data17,10)
data <- data17 %>%
  rowwise() %>%
  summarise(
    mean = mean(c_across(everything())),  
    sd = sd(c_across(everything()))     
  )
data<-as.data.frame(data)
data$K <- vector
data$pair <- c("GLUT-HK","HK-GPI", "GPI-PFK", "PFK-ALDO", "ALDO-TPI", "TPI-GAPDH", "GAPDH-PGK", "PGK-PGAM", "PGAM-ENO", "ENO-PK")
data$samp="KIRC"
summary<-rbind(summary,data)
data21 <- log(data21,10)
data <- data21 %>%
  rowwise() %>%
  summarise(
    mean = mean(c_across(everything())),  
    sd = sd(c_across(everything()))     
  )
data<-as.data.frame(data)
data$K <- vector
data$pair <- c("GLUT-HK","HK-GPI", "GPI-PFK", "PFK-ALDO", "ALDO-TPI", "TPI-GAPDH", "GAPDH-PGK", "PGK-PGAM", "PGAM-ENO", "ENO-PK")
data$samp="LUAD"
summary<-rbind(summary,data)
data25 <- log(data25,10)
data <- data25 %>%
  rowwise() %>%
  summarise(
    mean = mean(c_across(everything())),  
    sd = sd(c_across(everything()))     
  )
data<-as.data.frame(data)
data$K <- vector
data$pair <- c("GLUT-HK","HK-GPI", "GPI-PFK", "PFK-ALDO", "ALDO-TPI", "TPI-GAPDH", "GAPDH-PGK", "PGK-PGAM", "PGAM-ENO", "ENO-PK")
data$samp="LUSC"
summary<-rbind(summary,data)
data29 <- log(data29,10)
data <- data29 %>%
  rowwise() %>%
  summarise(
    mean = mean(c_across(everything())),  
    sd = sd(c_across(everything()))     
  )
data<-as.data.frame(data)
data$K <- vector
data$pair <- c("GLUT-HK","HK-GPI", "GPI-PFK", "PFK-ALDO", "ALDO-TPI", "TPI-GAPDH", "GAPDH-PGK", "PGK-PGAM", "PGAM-ENO", "ENO-PK")
data$samp="PRAD"
summary<-rbind(summary,data)
data33 <- log(data33,10)
data <- data33 %>%
  rowwise() %>%
  summarise(
    mean = mean(c_across(everything())),  
    sd = sd(c_across(everything()))     
  )
data<-as.data.frame(data)
data$K <- vector
data$pair <- c("GLUT-HK","HK-GPI", "GPI-PFK", "PFK-ALDO", "ALDO-TPI", "TPI-GAPDH", "GAPDH-PGK", "PGK-PGAM", "PGAM-ENO", "ENO-PK")
data$samp="THCA"
summary<-rbind(summary,data)


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

p<-ggplot(summary) + 
  geom_point(aes(y = K, x = mean, color = samp), size = 8, alpha = 0.5) +  
  geom_abline(
    data = line_df,
    aes(slope = 1, intercept = intercept, alpha = level),
    linetype = "dashed",
    color = "#4472C4"
  ) +
  scale_alpha_manual(
    name = "Difference between K and CAQ",
    values = c(
      "Equal" = 1,
      "10-fold difference" = 0.7,
      "100-fold difference" = 0.3
    )
  ) +
  guides(
    color = guide_legend(title = NULL, order = 1),
    alpha = guide_legend(title = NULL,override.aes = list(linetype = "dashed", color = "#4472C4", order = 2))
  ) +
  labs(
    title = expression(
      atop(
        "Comparison of K and CAQ for",
        "for glycolysis enzymes in human normal samples"
      )
    ),
    x = expression(log[10] * (CAQ)),  
    y = expression(log[10] * (K))
  ) +
  scale_color_manual(
    values = colour,
    labels = c("BRCA (n=112)", "COAD (n=41)", "HNSC (n=43)", "KIRC (n=71)", 
               "LUAD (n=58)", "LUSC (n=50)", "PRAD (n=51)", "THCA (n=58)"
               )  
  ) +
  theme(
    axis.text.x = element_text(size = 20,colour = "black"),  
    axis.title.x = element_text(size = 20,colour = "black"),  
    axis.text.y = element_text(size = 20,colour = "black"),  
    axis.title = element_text(size = 20,colour = "black"),  
    axis.ticks = element_line(colour = "black"),
    axis.line = element_line(colour = "black", linewidth = 0.8), 
    plot.title = element_text(size = 22, hjust = 0.5, vjust = 1), 
    plot.title.position = "plot",  
    panel.grid = element_blank(),  
    panel.background = element_rect(fill = "white", color = NA),  
    panel.border = element_rect(color = "black", fill = NA, size = 1),  
    legend.position = "right", 
    legend.background = element_rect(fill = "white", color = NA), 
    legend.text = element_text(size = 16)  
  ) +
  coord_cartesian(xlim = c(-6, 6), ylim = c(-6, 6)) +
  geom_text_repel(
    data = subset(summary, samp == "THCA"),
    aes(label = pair), size = 7,
    x = subset(summary, samp == "THCA")$mean,
    y = subset(summary, samp == "THCA")$K
  )
ggsave(paste0(dir,"/figure 4/figures/","figure 4b.pdf"), p, width = 800/96, height = 600/96, units = "in", device = "pdf")

