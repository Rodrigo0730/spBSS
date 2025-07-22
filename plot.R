#loads results from the analysis on puhti and performs a summary statistic of it

library(dplyr)
library(ggplot2)

#load results from supercomputer from e.g setting3 
results <- readRDS("setting3.rds")

results <- results[ !is.na(results$md) & !is.na(results$d), ]
summary <- results %>%
  group_by(model, d) %>%
  summarise(mean_md = mean(md), .groups = "drop")

ggplot(summary, aes(x = d, y = mean_md, color = model)) +
  geom_line() +
  geom_point() +
  labs(
    x = "d",
    y = "Average MDI value",
    color = "Model",
    title = ""
  ) +
  theme_minimal()
