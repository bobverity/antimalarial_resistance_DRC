# trends in COI and missingness
library(tidyverse)
library(rdhs)

figure_dir <- "analysis/manuscript_figures/"

# load our samples
samples <- readRDS("analysis/data/derived/samples_coi.rds")
samples$symptomatic <- TRUE
samples$symptomatic[samples$Country %in% c("DRC", "Zambia")] <- FALSE

source("analysis/facet_zoom.R")

# COI disstributions first
coi_dist <-
  ggplot(samples, aes(x = Country, y = coi)) + geom_violin(trim = FALSE) +
  stat_summary(
    fun.data = "mean_cl_boot",
    fun.args = list(B = 10000),
    geom = "pointrange",
    color = "red"
  ) +
  ylab("COI") + 
  scale_y_continuous(breaks = 0:8) + 
  theme(panel.grid.major.y = element_line(color = "grey"))

# Supp Figure 3
cowplot::save_plot(
  paste0(figure_dir, "/sup_fig3_v1.png"),
  coi_dist,
  base_height = 6,
  base_width = 9,
  dpi = 300
)
cowplot::save_plot(
  paste0(figure_dir, "/sup_fig3_v1.svg"),
  coi_dist,
  base_height = 6,
  base_width = 9,
  dpi = 300
)

cowplot::save_plot(
  paste0(figure_dir, "/sup_fig3_v2.png"),
  coi_dist + facet_zoom2(ylim = c(1.25, 3), zoom.size = 0.5),
  base_height = 6,
  base_width = 12,
  dpi = 300
)
cowplot::save_plot(
  paste0(figure_dir, "/sup_fig3_v2.svg"),
  coi_dist + facet_zoom2(ylim = c(1.25, 3), zoom.size = 0.5),
  base_height = 6,
  base_width = 12,
  dpi = 300
)

# and values of boostrap
boots <-
  group_by(samples, Country) %>% summarise(
    b = list(Hmisc::smean.cl.boot(coi, B = 10000)),
    mean = b[[1]][1],
    low = b[[1]][2],
    high = b[[1]][3]
  ) %>%
  select(Country, mean, low, high)

write.table(cbind(boots[, 1], round(boots[, 2:4], 3)),
            paste0(figure_dir, "/country_bootstrap_ci.csv"))

# Supp Figure 4
meta_admin <- readRDS("analysis/data/derived/admin_coi.rds")
meta_cluster <- readRDS("analysis/data/derived/cluster_coi.rds")

admin_plot <- ggplot(meta_admin, aes(micro, coi)) +
  geom_point(aes(size = n), alpha = 0.6) +
  geom_smooth(aes(weight = n), span = 3) +
  theme_bw() +
  xlab("Admin Microscopy Prevalence") +
  ylab("COI") + scale_size_continuous(name = "Sample Size")


cluster_plot <- ggplot(meta_cluster, aes(x = micro, y = coi, weight = n)) +
  geom_point(aes(size = n), alpha = 0.6) +
  geom_smooth(span = 1, se = TRUE) +
  theme_bw() + xlab("Cluster Microscopy Prevalence") +
  ylab("COI") + scale_size_continuous(name = "Sample Size")


cowplot::save_plot(
  paste0(figure_dir, "/sup_fig4_v1.png"),
  cowplot::plot_grid(admin_plot, cluster_plot, labels = "auto"),
  base_height = 6,
  base_width = 14,
  dpi = 300
)
cowplot::save_plot(
  paste0(figure_dir, "/sup_fig4_v1.svg"),
  cowplot::plot_grid(admin_plot, cluster_plot, labels = "auto"),
  base_height = 6,
  base_width = 14,
  dpi = 300
)

cowplot::save_plot(
  paste0(figure_dir, "/sup_fig4_v2.png"),
  cowplot::plot_grid(admin_plot, cluster_plot, labels = "auto"),
  base_height = 6,
  base_width = 14,
  dpi = 300
)
cowplot::save_plot(
  paste0(figure_dir, "/sup_fig4_v2.svg"),
  cowplot::plot_grid(admin_plot, cluster_plot, labels = "auto"),
  base_height = 6,
  base_width = 14,
  dpi = 300
)