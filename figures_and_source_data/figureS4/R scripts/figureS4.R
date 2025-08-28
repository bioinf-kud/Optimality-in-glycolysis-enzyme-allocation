dir<-"-path to folder-/figures_and_source_data"
head<-"/figureS4/data/"
data <- read.table(paste0(dir,head,"glyco_ratio.tsv"), header = TRUE, sep = "\t")
data$V2 <- factor(data$V2, levels = c("PYK", "ENO", "PGM", "PGK", "GAPD", "TPI", "FBA", "PFK", "PGI", "HEX1"))
ggplot(data, aes(x = Gene_name, y = mean, fill = V2)) +
  geom_bar(stat = "identity", width = 0.7) + 
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), 
                width = 0.2) +  
  theme_minimal() +
  labs(
    title = "Fraction of isoenzyme in TCGA samples",
    x = "Gene name",
    y = "Fraction of isoenzyme",
    fill = "BiGG reaction ID"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(
    values = c(
      "PYK" = "#B3CDE3", 
      "ENO" = "#A6D854", 
      "PGM" = "#8DA0CB", 
      "PGK" = "#E5C494", 
      "GAPD" = "#BEBADA", 
      "TPI" = "#B3B3B3", 
      "FBA" = "#F4CAE4", 
      "PFK" = "#FFF2AE", 
      "PGI" = "#A6D854", 
      "HEX1" = "#E5D8BD"
    )
  ) +
  facet_wrap(samp ~ V2, scales = "free_x", ncol = 10)

