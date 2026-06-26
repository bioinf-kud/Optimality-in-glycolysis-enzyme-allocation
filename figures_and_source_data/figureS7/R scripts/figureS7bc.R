library(ggplot2)
library(dplyr)

dir<-"-path to folder-/figures_and_source_data"
base <- "/figureS7/data/"

cancers <- c("TCGA-COAD","TCGA-BRCA","TCGA-HNSC","TCGA-KIRC",
             "TCGA-LUAD","TCGA-LUSC","TCGA-PRAD","TCGA-THCA")

suffix_map <- c(
  original = ".glycolysis.tsv",
  a_permutated = ".glycolysis.a_permuted.tsv",
  e_permutated = ".glycolysis.e_permuted.tsv",
  both_permutated = ".glycolysis.both_permuted.tsv"
)


vector <- c(5.522477554, 8.841552744, -0.6857715, 7.833065244, -0.772381588,
            -3.06580201, -6.172457509, 5.478969765, -1.774938001, 1.290864001)
vector <- vector / log(10)

read_and_prep <- function(cancer, suffix) {
  f <- paste0(dir, base, cancer, suffix)
  d <- read.table(f, sep = "\t", check.names = FALSE)
  

  if (cancer != "TCGA-COAD" && ncol(d) >= 10) {
    d <- d[, -12, drop = FALSE]
  }
  
  d <- t(d)
  colnames(d) <- make.unique(as.character(d[1, ]))
  d <- d[-1, , drop = FALSE]
  d <- as.data.frame(d, stringsAsFactors = FALSE)
  d
}


get_rmse <- function(data, vector, per_name, dat_name) {
  rmse_vals <- numeric(ncol(data))
  datanum <- numeric(ncol(data))
  
  for (i in seq_len(ncol(data))) {
    x <- log(as.numeric(data[[i]]),10)
    valid <- is.finite(x) & is.finite(vector)
    x <- x[valid]
    y <- vector[valid]
    if(length(x)<10){
      print(length(x))
      print(per_name)
      print(dat_name)
    }
    if (length(x) > 1) {
      datanum[i] <- length(x)
      rmse_vals[i] <- sqrt(mean((x - y)^2))
    } else {
      datanum[i] <- NA_real_
      rmse_vals[i] <- NA_real_
    }
  }
  
  data.frame(
    type_raw = names(data),
    value = rmse_vals,
    datanum = datanum,
    per = per_name,
    dat = dat_name,
    stringsAsFactors = FALSE
  )
}


all_list <- list()
k <- 1
for (ca in cancers) {
  for (nm in names(suffix_map)) {
    d <- read_and_prep(ca, suffix_map[[nm]])
    all_list[[k]] <- get_rmse(d, vector, per_name = nm, dat_name = ca)
    k <- k + 1
  }
}
dataf <- bind_rows(all_list)


dataf$type <- gsub("\\..*", "", dataf$type_raw)


dataf$per <- factor(
  dataf$per,
  levels = c("original", "a_permutated", "e_permutated", "both_permutated")
)


plot_rmse <- function(sample_type = c("Tumor", "Normal")) {
  sample_type <- match.arg(sample_type)
  

  cnt <- dataf %>%
    filter(type == sample_type, per == "original") %>%
    group_by(dat) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(new_label = paste0(dat, " (n=", count, ")"))
  label_map <- setNames(cnt$new_label, cnt$dat)
  
  data_sub <- dataf %>% filter(type == sample_type)
  title_line2 <- paste("glycolysis enzymes in human", tolower(sample_type) , "samples")
  ggplot(data_sub) +
    geom_boxplot(aes(y = value, x = per, fill = per, colour = per), alpha = 0.5) +
    scale_fill_manual(
      values = c("#824098","#70AD47","#ED7D31","#4472C4"),
      labels = c("original data", "a permuted", "[E] permuted", "both permuted")
    ) +
    scale_colour_manual(values = c("#824098","#70AD47","#ED7D31","#4472C4")) +
    labs(
      title = bquote(
        atop("RMSE between "* log[10]*K* " and " * log[10]*CAQ*" (EOI) for",.(title_line2))
      ),
      y = "RMSE",
      col = "method"
    ) +
    facet_wrap(~ dat, labeller = labeller(dat = label_map)) +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 20),
      axis.title = element_text(size = 20),
      plot.title = element_text(size = 22, hjust = 0.5, vjust = 1),
      plot.title.position = "plot",
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      legend.position = c(0.85, 0.17),
      legend.background = element_rect(fill = "white", color = NA),
      legend.text = element_text(size = 16),
      strip.text = element_text(size = 16)
    ) +
    guides(fill = guide_legend(title = NULL), colour = "none")
}


p_tumor <- plot_rmse("Tumor")
p_normal <- plot_rmse("Normal")

ggsave(paste0(dir,"/figureS7/figures/","figureS7b.pdf"), p_tumor, width = 800/96, height = 800/96, units = "in", device = "pdf")
ggsave(paste0(dir,"/figureS7/figures/","figureS7c.pdf"), p_normal, width = 800/96, height = 800/96, units = "in", device = "pdf")
