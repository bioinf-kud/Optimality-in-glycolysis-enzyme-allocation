library(dplyr)
library(broom)
dir<-"-path to folder-/figures_and_source_data"
head1<-"/figureS12/data/RMSE_list/TCGA-"
head2<-"/figureS12/data/pathway_matrix/TCGA-"
tail1<-".RMSE.tsv"
tail2<-"_KEGG_zscore.tsv"
corr1<-read.table(paste0(dir,head1,"BRCA",tail1),header=TRUE)
corr2<-read.table(paste0(dir,head1,"COAD",tail1),header=TRUE)
corr3<-read.table(paste0(dir,head1,"HNSC",tail1),header=TRUE)
corr4<-read.table(paste0(dir,head1,"KIRC",tail1),header=TRUE)
corr5<-read.table(paste0(dir,head1,"LUAD",tail1),header=TRUE)
corr6<-read.table(paste0(dir,head1,"LUSC",tail1),header=TRUE)
corr7<-read.table(paste0(dir,head1,"PRAD",tail1),header=TRUE)
corr8<-read.table(paste0(dir,head1,"THCA",tail1),header=TRUE)
kegg1<-read.table(paste0(dir,head2,"BRCA",tail2),header=TRUE,sep='\t')
kegg2<-read.table(paste0(dir,head2,"COAD",tail2),header=TRUE,sep='\t')
kegg3<-read.table(paste0(dir,head2,"HNSC",tail2),header=TRUE,sep='\t')
kegg4<-read.table(paste0(dir,head2,"KIRC",tail2),header=TRUE,sep='\t')
kegg5<-read.table(paste0(dir,head2,"LUAD",tail2),header=TRUE,sep='\t')
kegg6<-read.table(paste0(dir,head2,"LUSC",tail2),header=TRUE,sep='\t')
kegg7<-read.table(paste0(dir,head2,"PRAD",tail2),header=TRUE,sep='\t')
kegg8<-read.table(paste0(dir,head2,"THCA",tail2),header=TRUE,sep='\t')



data<-kegg1
data$corr<-corr1$value
data$sample <- ifelse(grepl("Tumor", rownames(data)), "Tumor",
                      ifelse(grepl("Normal", rownames(data)), "Normal", NA))
aggr<-data
aggr$type="BRCA"
aggrt<-aggr

data<-kegg2
data$corr<-corr2$value
data$sample <- ifelse(grepl("Tumor", rownames(data)), "Tumor",
                      ifelse(grepl("Normal", rownames(data)), "Normal", NA))
aggr<-data
aggr$type="COAD"
aggrt<-rbind(aggr,aggrt)

data<-kegg3
data$corr<-corr3$value
data$sample <- ifelse(grepl("Tumor", rownames(data)), "Tumor",
                      ifelse(grepl("Normal", rownames(data)), "Normal", NA))
aggr<-data
aggr$type="HNSC"
aggrt<-rbind(aggr,aggrt)

data<-kegg4
data$corr<-corr4$value
data$sample <- ifelse(grepl("Tumor", rownames(data)), "Tumor",
                      ifelse(grepl("Normal", rownames(data)), "Normal", NA))
aggr<-data
aggr$type="KIRC"
aggrt<-rbind(aggr,aggrt)

data<-kegg5
data$corr<-corr5$value
data$sample <- ifelse(grepl("Tumor", rownames(data)), "Tumor",
                      ifelse(grepl("Normal", rownames(data)), "Normal", NA))
aggr<-data
aggr$type="LUAD"
aggrt<-rbind(aggr,aggrt)

data<-kegg6
data$corr<-corr6$value
data$sample <- ifelse(grepl("Tumor", rownames(data)), "Tumor",
                      ifelse(grepl("Normal", rownames(data)), "Normal", NA))
aggr<-data
aggr$type="LUSC"
aggrt<-rbind(aggr,aggrt)

data<-kegg7
data$corr<-corr7$value
data$sample <- ifelse(grepl("Tumor", rownames(data)), "Tumor",
                      ifelse(grepl("Normal", rownames(data)), "Normal", NA))
aggr<-data
aggr$type="PRAD"
aggrt<-rbind(aggr,aggrt)

data<-kegg8
data$corr<-corr8$value
data$sample <- ifelse(grepl("Tumor", rownames(data)), "Tumor",
                      ifelse(grepl("Normal", rownames(data)), "Normal", NA))
aggr<-data
aggr$type="THCA"
aggrt<-rbind(aggr,aggrt)

group_colors <- c(
  "BRCA" = "#4472C4", 
  "COAD" = "#ED7D31", 
  "HNSC" = "#70AD47", 
  "KIRC" = "#824098",
  "LUAD" = "#C77EB5",
  "LUSC" = "#FFD166",
  "PRAD" = "#B89C6F",
  "THCA" = "#8BBDE3"
)

tumor_ribo <- aggrt[aggrt$sample == "Tumor", ]
r2_ribo <- tumor_ribo %>%
  filter(is.finite(corr), is.finite(hsa03010)) %>%
  group_by(type) %>%
  summarise(
    R2 = if (n() < 2) {
      NA_real_
    } else {
      broom::glance(stats::lm(corr ~ hsa03010, data = cur_data_all()))$r.squared
    },
    .groups = "drop"
  ) %>%
  mutate(
    R2_label = paste0(type, " (", round(R2, 3), ")")
  )
tumor_ribo <- tumor_ribo %>%
  left_join(r2_ribo, by = "type")
labels <- r2_ribo$R2_label
types <- r2_ribo$type
new_group_colors <- setNames(group_colors[types], labels)

p1<-ggplot(tumor_ribo, aes(x = hsa03010, y = corr, color = R2_label)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = new_group_colors, name = "Group (R²)")+
  labs(
    title = expression(
      atop(
        "Correlation between ribosome pathway activity",
        "and EOI of glycolysis in human tumor samples"
      )
    ),
    x = "Activity of ribosome pathway",
    y = "EOI of glycolysis"
  ) +
  theme(
    axis.text.x = element_text(hjust = 1, size = 20),  
    axis.text.y = element_text(size = 20),  
    axis.title = element_text(size = 20),  
    plot.title = element_text(size = 22, hjust = 0.5, vjust = 1), 
    plot.title.position = "plot",  
    panel.grid = element_blank(),  
    panel.background = element_rect(fill = "white", color = NA),  
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    legend.text = element_text(size = 16),        
    legend.title = element_text(size = 20) 
  )



tumor_oxphos <- aggrt[aggrt$sample == "Tumor", ]

r2_oxphos <- tumor_oxphos %>%
  filter(is.finite(corr), is.finite(hsa00190)) %>%
  group_by(type) %>%
  summarise(
    R2 = if (n() < 2) {
      NA_real_
    } else {
      broom::glance(stats::lm(corr ~ hsa00190, data = cur_data_all()))$r.squared
    },
    .groups = "drop"
  ) %>%
  mutate(
    R2_label = paste0(type, " (",round(R2, 3), ")")
  )
tumor_oxphos <- tumor_oxphos %>%
  left_join(r2_oxphos, by = "type")
labels <- r2_oxphos$R2_label
types <- r2_oxphos$type
new_group_colors <- setNames(group_colors[types], labels)

p2<-ggplot(tumor_oxphos, aes(x = hsa00190, y = corr, color = R2_label)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = new_group_colors, name = "Group (R²)") +
  labs(
    title = expression(
      atop(
        "Correlation between OXPHOS pathway activity",
        "and EOI of glycolysis in human tumor samples"
      )
    ),
    x = "Activity of OXPHOS pathway",
    y = "EOI of glycolysis"
  ) +
  theme(
    axis.text.x = element_text(hjust = 1, size = 20),  
    axis.text.y = element_text(size = 20),  
    axis.title = element_text(size = 20),  
    plot.title = element_text(size = 22, hjust = 0.5, vjust = 1), 
    plot.title.position = "plot",  
    panel.grid = element_blank(),  
    panel.background = element_rect(fill = "white", color = NA),  
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    legend.text = element_text(size = 16),        
    legend.title = element_text(size = 20) 
  )
ggsave(paste0(dir,"/figureS12/figures/","figureS12b.pdf"), p1, width = 800/96, height = 600/96, units = "in", device = "pdf")
ggsave(paste0(dir,"/figureS12/figures/","figureS12d.pdf"), p2, width = 800/96, height = 600/96, units = "in", device = "pdf")






