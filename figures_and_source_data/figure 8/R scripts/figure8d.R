library(dplyr)
library(ggplot2)
library(ggrepel) 
library(ggbreak)
dir<-"-path to folder-/figures_and_source_data"
file <- paste0(dir, "/figure 8/data/corr_P.tsv")
data <- read.table(file, header = TRUE, sep = "\t")


a <- 9
b <- 20
gap <- b - a
visual_gap <- 0.6 

data2 <- data %>%
  mutate(
    y_raw = -log10(p),
    y_plot = case_when(
      y_raw <= a ~ y_raw,
      y_raw >= b ~ y_raw - gap + visual_gap,
      TRUE ~ NA_real_
    )
  )

y_breaks_raw <- c(0, 3, 6, 9, 20, 21)
y_breaks_raw <- y_breaks_raw[y_breaks_raw <= max(data2$y_raw, na.rm = TRUE)]
y_breaks_plot <- ifelse(
  y_breaks_raw <= a,
  y_breaks_raw,
  y_breaks_raw - gap + visual_gap
)

sig_raw <- -log10(0.05)
sig_plot <- if (sig_raw <= a) {
  sig_raw
} else if (sig_raw >= b) {
  sig_raw - gap + visual_gap
} else {
  NA_real_
}

p <- ggplot(subset(data2, !is.na(y_plot)),
            aes(x = rho, y = y_plot, color = species, label = type)) +
  geom_point(alpha = 0.5, size = 5) +
  geom_text_repel(
    data = subset(data2, p < 0.05 & !is.na(y_plot)),
    size = 5,
    max.overlaps = 30,
    box.padding = 0.5,
    segment.color = "grey50"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#824098", size = 1, alpha = 0.7) +
  geom_hline(yintercept = sig_plot, linetype = "dashed", color = "#824098", size = 1, alpha = 0.7) +
  scale_y_continuous(
    breaks = y_breaks_plot,
    labels = y_breaks_raw
  ) +
  scale_color_manual(
    values = c(
      "H.sapiens" = "#4472C4",
      "E.coli" = "#ED7D31",
      "S.cerevisiae" = "#70AD47"
    ),
    labels = c(
      "H.sapiens" = expression(italic("H.sapiens")),
      "E.coli" = expression(italic("E.coli")),
      "S.cerevisiae" = expression(italic("S.cerevisiae"))
    ),
    name = "species"
  ) +
  labs(
    x = expression("Spearman's " ~ rho),
    y = expression(-log[10] * "(p-value)"),
    color = "Experimental condition",
    title = "Correlation between glycolysis enzyme\n fraction and EOI of glycolysis"
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
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18),
    legend.position = c(0.20, 0.60)
  )+
  coord_cartesian(xlim = c(-1, 1), clip = "off")

p <- p +
  annotate(
    "segment",
    x = -1.15, xend = -1.05,
    y = 9.1, yend = 9.3,
    linewidth = 1,
    color = "black"
  ) +
  annotate(
    "segment",
    x = -1.15, xend = -1.05,
    y = 9.3, yend = 9.5,
    linewidth = 1,
    color = "black"
  )
p <- p +
  annotate(
    "segment",
    x =-0.05, xend = 0.05,
    y = 9.1, yend = 9.3,
    linewidth = 1,
    color = "#824098", size = 1, alpha = 0.7
  ) +
  annotate(
    "segment",
    x =-0.05, xend = 0.05,
    y = 9.3, yend = 9.5,
    linewidth = 1,
    color = "#824098", size = 1, alpha = 0.7
  )
ggsave(
  paste0(dir, "/figure 5/figures/", "figure 5d.pdf"),
  p, width = 650/96, height = 700/96, units = "in", device = "pdf"
)

ggsave(paste0(dir,"/figure 8/figures/","figure 8d.pdf"), p, width = 650/96, height = 700/96, units = "in", device = "pdf")

