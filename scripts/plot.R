#loads results from the analysis on puhti and performs a summary statistic of it

library(dplyr)
library(ggplot2)

#load results from supercomputer from e.g setting3
#better to separate the results per value of d, that way
#if code in puhti fails for memory issues we can at least
#recover the results from which it last failed.

#results are in /results/raw
results <- readRDS("results/raw/...")

results <- results[ !is.na(results$md) & !is.na(results$d), ]
summary <- results %>%
  group_by(model, d) %>%
  summarise(mean_md = mean(md), .groups = "drop")

plot <- ggplot(summary, aes(x = d, y = mean_md, color = model)) +
  geom_line() +
  geom_point() +
  labs(
    x = "d",
    y = "Average MDI value",
    color = "Model",
    title = ""
  ) +
  theme_minimal()

filename <- "/results/plots/setting_3"
ggsave(filename, plot = plot, width = 6, height = 4, dpi = 300)



