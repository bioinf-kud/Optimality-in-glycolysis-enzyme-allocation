library(ggplot2)
library(dplyr)
library(tidyr)

dir<-"-path to folder-/figures_and_source_data"
head <- "/figure 4/data/calculated/"
tail1 <- ".glycolysis.tsv"

vector <- c(5.522477554, 8.841552744, -0.6857715, 7.833065244, -0.772381588,
            -3.06580201, -6.172457509, 5.478969765, -1.774938001, 1.290864001)
vector <- vector / log(10)
read_and_process <- function(filename, is_brca = FALSE) {
  data <- read.table(filename, sep = "\t", check.names = FALSE)
  if (is_brca) {
    data <- data[, -12, drop = FALSE]
  }
  data <- t(data)
  colnames(data) <- make.unique(as.character(data[1, ]))
  data <- data[-1, , drop = FALSE]
  as.data.frame(data, stringsAsFactors = FALSE)
}

datasets <- list(
  "TCGA-COAD" = read_and_process(paste0(dir, head, "TCGA-COAD", tail1)),
  "TCGA-BRCA" = read_and_process(paste0(dir, head, "TCGA-BRCA", tail1), is_brca = TRUE),
  "TCGA-HNSC" = read_and_process(paste0(dir, head, "TCGA-HNSC", tail1), is_brca = TRUE),
  "TCGA-KIRC" = read_and_process(paste0(dir, head, "TCGA-KIRC", tail1), is_brca = TRUE),
  "TCGA-LUAD" = read_and_process(paste0(dir, head, "TCGA-LUAD", tail1), is_brca = TRUE),
  "TCGA-LUSC" = read_and_process(paste0(dir, head, "TCGA-LUSC", tail1), is_brca = TRUE),
  "TCGA-PRAD" = read_and_process(paste0(dir, head, "TCGA-PRAD", tail1), is_brca = TRUE),
  "TCGA-THCA" = read_and_process(paste0(dir, head, "TCGA-THCA", tail1), is_brca = TRUE)
)

get_rmse <- function(data, vector) {
  rmse_vals <- numeric(ncol(data))
  types <- gsub("\\..*", "", names(data))
  
  for (i in seq_len(ncol(data))) {
    current_column <- log(as.numeric(data[[i]]),10)
    valid <- is.finite(current_column) & is.finite(vector)
    x <- current_column[valid]
    y <- vector[valid]
    
    if (length(x) > 1) {
      rmse_vals[i] <- sqrt(mean((x - y)^2))
    } else {
      rmse_vals[i] <- NA_real_
    }
  }
  
  data.frame(
    RMSE = rmse_vals,
    Group = types,
    stringsAsFactors = FALSE
  )
}

tcga_rmse <- lapply(names(datasets), function(dat_name) {
  df <- get_rmse(datasets[[dat_name]], vector)
  df$Dataset <- dat_name
  df
}) %>% bind_rows()

tcga_rmse <- tcga_rmse %>% filter(Group %in% c("Tumor", "Normal"))

prot_head <- "/figure 4/data/calculated/CCLE_proteome.glycolysis"
prot_tail <- ".tsv"

read_proteome <- function(filename) {
  data <- read.table(filename, sep = "\t", check.names = FALSE)
  data<-subset(data,select=-V12)
  rn <- make.unique(as.character(data$V1))
  mat <- as.matrix(data[, -1, drop = FALSE])
  rownames(mat) <- rn
  as.data.frame(t(mat), stringsAsFactors = FALSE)
}

proteome_data <- read_proteome(paste0(dir, prot_head, prot_tail))
proteome_rmse <- get_rmse(proteome_data, vector)
proteome_rmse$Group <- "Proteome"

merged <- bind_rows(
  tcga_rmse %>% select(Group, RMSE),
  proteome_rmse %>% select(Group, RMSE)
)

merged$Group <- factor(
  merged$Group,
  levels = c("Tumor", "Normal", "Proteome"),
  labels = c("TCGA tumor", "TCGA normal", "CCLE cell lines")
)

p<-ggplot(merged, aes(x = Group, y = RMSE)) +
  geom_violin(
    fill = "#AFCBE3",
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
  labs(
    y = "RMSE",
    title = expression(
        atop("RMSE between "* log[10]*K* " and " * log[10]*CAQ*" (EOI)","for glycolysis enzymes in human")
    )
  ) +
  xlab("")
ggsave(paste0(dir,"/figure 4/figures/","figure 4e.pdf"), p, width = 650/96, height = 700/96, units = "in", device = "pdf")

