library(ggplot2)
library(dplyr)
library(tidyr)


dir<-"-path to folder-/figures_and_source_data"
head <- "/figureS5/data/calculated/"
tail1 <- ".glycolysis.tsv"
vector <- c(5.522477554, 8.841552744, -0.6857715, 7.833065244, -0.772381588,
            -3.06580201, -6.172457509, 5.478969765, -1.774938001, 1.290864001)
vector <- vector / log(10)
read_and_process <- function(filename, is_brca=FALSE) {
  data <- read.table(filename, sep="\t")
  if (is_brca) {
    data <- data[, -12]
  }
  data <- t(data)
  colnames(data) <- data[1, ]
  data <- data[-1, ]
  as.data.frame(data)
}

datasets <- list(
  "TCGA-COAD" = read_and_process(paste0(dir, head, "TCGA-COAD", tail1)),
  "TCGA-BRCA" = read_and_process(paste0(dir, head, "TCGA-BRCA", tail1), is_brca=TRUE),
  "TCGA-HNSC" = read_and_process(paste0(dir, head, "TCGA-HNSC", tail1), is_brca=TRUE),
  "TCGA-KIRC" = read_and_process(paste0(dir, head, "TCGA-KIRC", tail1), is_brca=TRUE),
  "TCGA-LUAD" = read_and_process(paste0(dir, head, "TCGA-LUAD", tail1), is_brca=TRUE),
  "TCGA-LUSC" = read_and_process(paste0(dir, head, "TCGA-LUSC", tail1), is_brca=TRUE),
  "TCGA-PRAD" = read_and_process(paste0(dir, head, "TCGA-PRAD", tail1), is_brca=TRUE),
  "TCGA-THCA" = read_and_process(paste0(dir, head, "TCGA-THCA", tail1), is_brca=TRUE)
)


get_pearson <- function(data, vector) {
  correlations <- numeric(ncol(data))
  types <- gsub("\\..*", "", names(data))
  for (i in 1:ncol(data)) {
    current_column <- as.numeric(data[[i]])
    current_column <- log(current_column)
    valid_indices <- which(current_column != -Inf & current_column != Inf)
    filtered_column <- current_column[valid_indices]
    filtered_vector <- vector[valid_indices]
    if (length(valid_indices) > 1) {
      correlations[i] <- cor(filtered_column, filtered_vector, method = "pearson")
    } else {
      correlations[i] <- NA
    }
  }
  data.frame(
    pearsonCorrelation = correlations,
    Group = types
  )
}

tcga_corr <- lapply(names(datasets), function(dat_name) {
  df <- get_pearson(datasets[[dat_name]], vector)
  df$Dataset <- dat_name
  df
}) %>% bind_rows()

tcga_corr <- tcga_corr %>% filter(Group %in% c("Tumor", "Normal"))

prot_head <- "/figure 4/data/calculated/CCLE_proteome.glycolysis"
prot_tail <- ".tsv"
read_proteome <- function(filename) {
  data <- read.table(filename, sep="\t")
  row.names(data) <- data$V1
  data <- subset(data, select = -V1)
  t(data) %>% as.data.frame()
}

proteome_data <- read_proteome(paste0(dir, prot_head, prot_tail))
proteome_corr <- get_pearson(proteome_data, vector)
proteome_corr$Group <- "Proteome"


merged <- dplyr::bind_rows(
  dplyr::select(tcga_corr, Group, pearsonCorrelation),
  dplyr::select(proteome_corr, Group, pearsonCorrelation)
)
merged$Group <- factor(merged$Group,
                       levels = c("Tumor", "Normal", "Proteome"),
                       labels = c("TCGA tumor", "TCGA normal", "CCLE cell lines")
)

p<-ggplot(merged, aes(x = Group, y = pearsonCorrelation)) +
  geom_violin(
    fill = '#AFCBE3',
    color = "#4472C4",
    alpha = 0.6,
    draw_quantiles = c(0.5)
  ) +
  theme(
    axis.text.x = element_text(angle = 10, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 22, hjust = 0.5, vjust = 1),
    plot.title.position = "plot",
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  labs(y = expression("Pearson's r"),
       title = expression(
         atop(
           "Correlation between K and CAQ (EOI)",
           "for glycolysis enzymes in human"
         )
       )
       ) +
  xlab("") +
  ylim(-1, 1)   
ggsave(paste0(dir,"/figureS5/figures/","figureS5a.pdf"), p, width = 650/96, height = 700/96, units = "in", device = "pdf")
