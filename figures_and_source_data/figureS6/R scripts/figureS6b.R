library(ggplot2)
library(ggrepel)

dir<-"-path to folder-/figures_and_source_data"
head <- "/figureS6/data/"

data <- read.table(paste0(dir, head, "sce.glycolysis.tsv"), sep = "\t", check.names = FALSE)
data1 <- read.table(paste0(dir, head, "sce.glycolysis.no_substrate.tsv"), sep = "\t", check.names = FALSE)

vector <- c(
  5.522477554, 8.111677106, -0.6857715, 7.103189605, -0.29214456,
  -3.065802002, -5.00598449, 6.209460215, -1.774938001, 1.290864001
)
vector <- vector / log(10)
row.names(data) <- make.unique(as.character(data$V1))
row.names(data1) <- make.unique(as.character(data1$V1))

data <- subset(data, select = -V1)
data <- subset(data, select = -V10)
data1 <- subset(data1, select = -V1)
data1 <- subset(data1, select = -V10)

data <- as.data.frame(t(data), stringsAsFactors = FALSE)
data1 <- as.data.frame(t(data1), stringsAsFactors = FALSE)

log_data <- as.data.frame(lapply(data, function(x) log(as.numeric(x),10)))
log_data1 <- as.data.frame(lapply(data1, function(x) log(as.numeric(x,10))))

rmse <- numeric(ncol(log_data))
rmse1 <- numeric(ncol(log_data1))

for (i in seq_len(ncol(log_data))) {
  x <- as.numeric(log_data[[i]])
  valid <- is.finite(x) & is.finite(vector)
  rmse[i] <- if (sum(valid) > 1) sqrt(mean((x[valid] - vector[valid])^2)) else NA_real_
}

for (i in seq_len(ncol(log_data1))) {
  x <- as.numeric(log_data1[[i]])
  valid <- is.finite(x) & is.finite(vector)
  rmse1[i] <- if (sum(valid) > 1) sqrt(mean((x[valid] - vector[valid])^2)) else NA_real_
}

rmse_df <- data.frame(value = rmse, method = "with_substrate")
row.names(rmse_df) <- names(data)

rmse_df1 <- data.frame(value = rmse1, method = "without_substrate")
row.names(rmse_df1) <- names(data1)

Amean <- mean(rmse_df$value, na.rm = TRUE)
Bmean <- mean(rmse_df1$value, na.rm = TRUE)

plot_df <- rbind(rmse_df, rmse_df1)

p<-ggplot(plot_df) +
  ggtitle("RMSE distribution (n=16)")+
  geom_histogram(aes(x = value, fill = method, colour = method),
                 binwidth = 0.05, position = "identity", alpha = 0.5) +
  scale_fill_manual(
    values = c("#4472C4", "#ED7D31"),
    labels = c(
      paste0("Adjusted a (mean RMSE: ", sprintf("%.2f", Amean), ")"),
      paste0("Original a (mean RMSE: ", sprintf("%.2f", Bmean), ")")
    )
  ) +
  scale_colour_manual(values = c("#4472C4", "#ED7D31")) +
  labs(
    title = expression(
      atop("RMSE between "* log[10]*K* " and " * log[10]*CAQ*" (EOI) for","glycolysis enzymes in 16 " * italic("S.cerevisiae") * " samples")
    ),
    x = "RMSE",
    y = "Number of samples",
    col = "method"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 22, hjust = 0.5, vjust = 1),
    plot.title.position = "plot",
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    legend.position = c(0.69, 0.94),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.text = element_text(size = 16)
  ) +
  guides(fill = guide_legend(title = NULL), colour = "none")
ggsave(paste0(dir,"/figureS6/figures/","figureS6b.pdf"), p, width = 650/96, height = 700/96, units = "in", device = "pdf")
