library(dplyr)
library(ggplot2)

dir<-"-path to folder-/figures_and_source_data"
types <- c("BRCA", "COAD", "HNSC", "KIRC", "LUAD", "LUSC", "PRAD", "THCA")

all_data <- list()
for (type in types) {
  file <- paste0(dir, "/figureS10/data/glycolysis_tpm_sum/TCGA-", type, "_glyco_sum_RMSE_Tumor.tsv")
  df <- read.table(file, header = TRUE, sep = "\t")
  df$type <- type
  all_data[[type]] <- df
}

aggrt <- do.call(rbind, all_data)
aggrt$Glyco_Genes_Sum <- log(aggrt$Glyco_Genes_Sum)

corr_stats <- aggrt %>%
  group_by(type) %>%
  summarise(
    rho = cor(Glyco_Genes_Sum, RMSE, method = "spearman", use = "complete.obs"),
    p = cor.test(Glyco_Genes_Sum, RMSE, method = "spearman")$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    Rho_label = paste0(type, " (", round(rho, 2), ", ", formatC(p, format = "e", digits = 2), ")")
  )

aggrt <- aggrt %>% left_join(corr_stats %>% select(type, Rho_label), by = "type")

aggr_df <- aggrt %>%
  group_by(type, RMSE) %>%
  summarise(
    Glyco_Genes_Sum = mean(Glyco_Genes_Sum, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(corr_stats %>% select(type, Rho_label), by = "type")

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
new_group_colors <- setNames(group_colors[corr_stats$type], corr_stats$Rho_label)

p<-ggplot(aggrt, aes(x = Glyco_Genes_Sum, y = RMSE, color = Rho_label)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) + 
  scale_color_manual(values = new_group_colors, name = expression("Cancer type (" ~ rho ~ ", p-value)")) +
  facet_wrap(~ type, scales = "free") +
  labs(
    x = "Glycolysis enzyme expression (log(TPM))",
    y = "EOI of glycolysis",
    color = "Cancer type",
    title = "Expression of glycolysis enzyme and EOI\nof glycolysis in human tumor samples"
  ) +
  theme(
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(hjust = 1, size = 16),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 22, hjust = 0.5, vjust = 1),
    plot.title.position = "plot",
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    legend.position = c(0.84, 0.17),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 20),
    strip.text = element_text(size = 16)
  )
ggsave(paste0(dir,"/figureS10/figures/","figureS10b.pdf"), p, width = 1000/96, height = 1000/96, units = "in", device = "pdf")

