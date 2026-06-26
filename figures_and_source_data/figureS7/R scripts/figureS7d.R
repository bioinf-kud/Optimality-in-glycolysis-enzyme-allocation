library(ggplot2)
library(ggrepel)

dir<-"-path to folder-/figures_and_source_data"
head <- "/figureS7/data/eco.glycolysis"
tail <- ".tsv"

data1 <- read.table(paste0(dir, head, tail), sep = "\t")
data2 <- read.table(paste0(dir, head, ".a_permuted", tail), sep = "\t")
data3 <- read.table(paste0(dir, head, ".e_permuted", tail), sep = "\t")
data4 <- read.table(paste0(dir, head, ".both_permuted", tail), sep = "\t")

prep <- function(d) {
  d <- subset(d, select = -V1)
  as.data.frame(t(d))
}
prep1 <- function(d) {
  d <- subset(d, select = -V1)
  d <- subset(d, select = -V11)
  as.data.frame(t(d))
}

data1 <- prep1(data1)
data2 <- prep1(data2)

data3 <- prep(data3)
data4 <- prep(data4)

vector <- c(
  13.34312772, -0.077827744, 8.924192273, -0.29214456, -2.006169471,
  -3.01912952, 3.362797029, -1.682640295, 1.171297752
)
vector <- vector / log(10)

calc_rmse_df <- function(data, vector, label) {
  vals <- numeric(ncol(data))
  datanum <- numeric(ncol(data))
  
  for (i in seq_len(ncol(data))) {
    current_column <- log(as.numeric(data[[i]]),10)
    valid <- is.finite(current_column) & is.finite(vector)
    
    x <- current_column[valid]
    y <- vector[valid]
    
    if (length(x) > 1) {
      datanum[i] <- length(x)
      vals[i] <- sqrt(mean((x - y)^2))
    } else {
      vals[i] <- NA_real_
    }
  }
  
  out <- data.frame(value = vals, datanum = datanum, per = label)
  out
}

d1 <- calc_rmse_df(data1, vector, "original data")
d2 <- calc_rmse_df(data2, vector, "a permuted")
d3 <- calc_rmse_df(data3, vector, "[E] permuted")
d4 <- calc_rmse_df(data4, vector, "both permuted")
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
      atop("RMSE between "* log[10]*K* " and " * log[10]*CAQ*" (EOI) for","glycolysis enzymes in " * italic("E.coli")*" samples")
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
ggsave(paste0(dir,"/figureS7/figures/","figureS7d.pdf"), p, width = 650/96, height = 700/96, units = "in", device = "pdf")
