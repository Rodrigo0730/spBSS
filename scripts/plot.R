# This is the file which produced the figures found in the paper. 
# Creates each plot from each setting independently and then arranges
# everything together to obtain the final figure.

library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(cowplot)


load_means <- function(setting, ds) {
  lapply(ds, function(d) {
    res <- readRDS(sprintf("results/raw/setting_%s/res_%d.rds", setting, d))
    m   <- colMeans(res, na.rm = TRUE)
    data.frame(
      setting = setting,
      d       = d,
      n       = d^2,
      model   = c("spFOBI", "spJADE", "SBSS", "FOBI", "JADE"),
      mean    = as.numeric(m),
      stringsAsFactors = FALSE
    )
  }) %>%
    bind_rows()
}

ds <- c(5, 10, 20, 40, 60, 80)
all_means <- bind_rows(
  load_means(1, ds),
  load_means(2, ds),
  load_means(3, ds),
  load_means(4, ds)
)

palette <- brewer.pal(5, "Set2")

theme_journal <- function() {
  theme_classic(base_size = 12, base_family = "serif") +
    theme(
      panel.border       = element_rect(color = "black", fill = NA, size = 0.5),
      panel.grid.major.y = element_line(color = "grey85", size = 0.4),
      axis.text          = element_text(color = "black", size = 8),
      axis.title         = element_text(size = 11),
      axis.line          = element_line(linewidth = 0.1),
      plot.title         = element_text(size = 10, hjust = 0.5, vjust=0, face = "bold"),
      legend.position    = "bottom",
      legend.title       = element_blank(),
      legend.key         = element_blank(),
      legend.background  = element_rect(fill = alpha("white", 0.7), color = NA),
      strip.background   = element_blank(),
      strip.text         = element_text(face = "bold")
    )
}

plot_setting <- function(df, setting) {
  ggplot(df, aes(
    x     = factor(d),
    y     = mean,
    color = model,
    group = model
  )) +
    geom_line(size = 0.7, na.rm = TRUE) +
    geom_point(size = 0.5, na.rm = TRUE) +
    scale_x_discrete(
      name   = "n",
      labels = function(x) as.numeric(x)^2,
      expand = c(0.03, 0.03)
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      expand = c(0.03, 0.03)
    ) +
    scale_color_manual(values = palette) +
    labs(
      title = paste("Setting", setting),
      y     = ""
    ) +
    theme_journal()
}

placeholder_plot <- function(setting) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = paste("Setting", setting, "\n(no data yet)"),
             size = 4, fontface = "italic") +
    xlim(0, 1) + ylim(0, 1) +
    theme_void(base_size = 12, base_family = "serif") +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      plot.title   = element_text(size = 12, hjust = 0.5, face = "bold")
    ) +
    labs(title = paste("Setting", setting))
}
# p1 <- placeholder_plot(1)
p1 <- plot_setting(all_means %>% filter(setting == 1), 1) +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
p2 <- plot_setting(all_means %>% filter(setting == 2), 2) +
      theme(axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank()
            )
      
p3 <- plot_setting(all_means %>% filter(setting == 3), 3) + 
      theme(axis.title.x = element_blank())

p4 <- plot_setting(all_means %>% filter(setting == 4), 4) +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            axis.title.x = element_blank())

xlabel <- ggplot() + 
  theme_void() + 
  annotate("text",
           x = 0.5,
           y = 0.5,
           label = "Number of spatial locations",
           size = 3,
           family = "serif") + 
  theme(plot.margin = margin())

ylabel <- ggplot() +
  theme_void() +
  annotate("text",
           x = 0.5, y = 0.5,
           label = "Average of Minimum Distance Index",
           angle = 90,
           size = 3,
           family = "serif") + 
  theme(plot.margin = margin())


combined <- ((p1 + p2) / (p3 + p4) / xlabel) +
  plot_annotation(title = NULL,
                  caption = NULL,
                  tag_levels = NULL) +
  plot_layout(heights = c(1, 1, 0.2),
              guides = "collect") &
  theme(legend.position = "bottom",
        plot.margin = margin())

res <- (ylabel | combined) + 
  plot_layout(widths = c(0.02, 1))

print(res)



