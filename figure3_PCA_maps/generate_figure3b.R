
# generate_figure3b.R
#
# Author: Bob Verity
# Date: 2019-12-17
#
# Purpose:
# Produce plots of PrevMap model predictions from drug resistance frequencies.
#
# ------------------------------------------------------------------

# load packages
library(ggplot2)
library(rworldmap)
library(viridisLite)
library(gridExtra)

# ------------------------------------------------------------------

# define plotting parameters
col_country <- grey(0.95)
col_country_border <- grey(0.5)
size_country_border <- 0.5
#col_sea <- grey(1.0)
shape_resolution <- "coarse"
#col_limits <- c(-4,4)
#point_size <- 2.5
#stroke_col <- grey(0.5)
#stroke_size <- 0.25
col_vec <- viridisLite::plasma(100)

# read in predictions
pred <- readRDS("generated_data/dr_PrevMap_predictions.rds")

# load country shapefiles, leaving gap for DRC
world_map <- getMap(resolution = shape_resolution)
sub_world_map <- subset(world_map, ISO3 != "COD")

# create basic map plot
plot_base <- ggplot() + theme_bw() +
  theme(panel.background = element_rect(fill = col_sea),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(angle = 0, hjust = 0, size = 12))

# create list of plot objects
plot_list <- list()
for (i in 1:length(fit)) {
  
  plot_list[[i]] <- plot_base +
    geom_raster(aes(x = x, y = y, fill = layer), data = pred[[i]]$df_rast) +
    geom_contour(aes(x = x, y = y, z = layer), size = 0.5, color = "black", data = pred[[i]]$df_rast_smooth) +
    geom_polygon(aes(long, lat, group = group),
                 size = size_country_border, color = col_country_border,
                 fill = col_country, data = sub_world_map) +
    coord_cartesian(xlim = c(11, 33), ylim = c(-13,6)) +
    scale_fill_gradientn(colours = col_vec, name = "prevalence") +
    xlab("longitude") + ylab("latitude")
}

# finalise plots
plot1 <- plot_list[[1]] + ggtitle(expression(paste("e) ", italic("pfcrt"), " K76T prevalence")))
plot2 <- plot_list[[2]] + ggtitle(expression(paste("f) ", italic("dhps"), " K540E prevalence")))

# produce combined plot
plot_combined <- gridExtra::grid.arrange(plot1, plot2, ncol = 1)

# save to file
ggsave("figure3_PCA_maps/figure3b_DR_PrevMaps.pdf", plot = plot_combined,
       device = "pdf", width = 4.5, height = 7)
ggsave("figure3_PCA_maps/figure3b_DR_PrevMaps.png", plot = plot_combined,
       device = "png", width = 4.5, height = 7, dpi = 100)

