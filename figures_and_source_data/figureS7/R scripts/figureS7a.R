library(ggplot2)
library(ggrepel)

dir<-"-path to folder-/figures_and_source_data"
head <- "/figureS7/data/CCLE_proteome.glycolysis"
tail <- ".tsv"

data1 <- read.table(paste0(dir, head, tail), sep = "\t")
data2 <- read.table(paste0(dir, head, ".a_permuted", tail), sep = "\t")
data3 <- read.table(paste0(dir, head, ".e_permuted", tail), sep = "\t")
data4 <- read.table(paste0(dir, head, ".both_permuted", tail), sep = "\t")

prep <- function(d) {
  row.names(d) <- d$V1
  d <- subset(d, select = -V1)
  d <- as.data.frame(t(d))
  d
}

prep1 <- function(d) {
  row.names(d) <- d$V1
  d <- subset(d, select = -V1)
  d <- subset(d, select = -V12)
  d <- as.data.frame(t(d))
  d
}
data1 <- prep1(data1)
data2 <- prep1(data2)
data3 <- prep(data3)
data4 <- prep(data4)

vector <- c(5.522477554, 8.841552744, -0.6857715, 7.833065244, -0.772381588,
            -3.06580201, -6.172457509, 5.478969765, -1.774938001, 1.290864001)
vector <- vector / log(10)
get_rmse <- function(data, vector, per_name) {
  rmse_vals <- numeric(ncol(data))
  datanum <- numeric(ncol(data))
  
  for (i in seq_len(ncol(data))) {
    current_column <- as.numeric(data[[i]])
    current_column <- log(current_column,10)
    
    valid <- is.finite(current_column) & is.finite(vector)
    x <- current_column[valid]
    y <- vector[valid]
    
    if (length(x) > 1) {
      datanum[i] <- length(x)
      rmse_vals[i] <- sqrt(mean((x - y)^2))
    } else {
      rmse_vals[i] <- NA_real_
    }
  }
  
  out <- data.frame(
    value = rmse_vals,
    datanum = datanum,
    per = per_name
  )
  row.names(out) <- names(data)
  out
}

d1 <- get_rmse(data1, vector, "original data")
d2 <- get_rmse(data2, vector, "a permuted")
d3 <- get_rmse(data3, vector, "[E] permuted")
d4 <- get_rmse(data4, vector, "both permuted")

dataf <- rbind(d1, d2, d3, d4)
dataf$per <- factor(dataf$per, levels = c("original data", "a permuted", "[E] permuted", "both permuted"))

p<-ggplot(dataf) +
  geom_boxplot(aes(y = value, x = per, fill = per, colour = per), alpha = 0.5) +
  geom_jitter(aes(y = value, x = per, colour = per), size = 1, alpha = 0.1) +
  scale_fill_manual(values = c(
    "both permuted" = "#4472C4",
    "[E] permuted" = "#ED7D31",
    "a permuted" = "#70AD47",
    "original data" = "#824098"
  )) +
  scale_colour_manual(values = c(
    "both permuted" = "#4472C4",
    "[E] permuted" = "#ED7D31",
    "a permuted" = "#70AD47",
    "original data" = "#824098"
  )) +
  labs(
    title = expression(
      atop("RMSE between "* log[10]*K* " and " * log[10]*CAQ*" (EOI) for","glycolysis enzymes in 378 human cell lines")
    ),
    y = "RMSE",
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
    legend.position = "none"
  ) +
  guides(fill = guide_legend(title = NULL), colour = "none")
ggsave(paste0(dir,"/figureS7/figures/","figureS7a.pdf"), p, width = 650/96, height = 700/96, units = "in", device = "pdf")

