
# generate_figure3.R
#
# Author: Bob Verity
# Date: 2019-12-17
#
# Purpose:
# Produce maps of first few principal components.
#
# ------------------------------------------------------------------

# NOTE - uncomment these lines to install packages as needed
#devtools::install_github("mrc-ide/MIPAnalyzer", ref = "version1.0.0")
#devtools::install_github("bobverity/bobfunctions2", ref = "version1.0.0")

# load packages
library(devtools)
library(MIPanalyzer)
library(tidyr)
library(bobfunctions2)
library(rworldmap)

# ------------------------------------------------------------------

# read in raw data
dat <- readRDS("source_data/biallelic_distances.rds")

# PCA on WSAFs
wsaf_impute <- get_wsaf(dat, impute = TRUE, FUN = mean)
pca <- pca_wsaf(wsaf_impute)

# find which points belong in which spatial clusters
dat$samples$lonlat <- paste(dat$samples$long, dat$samples$lat, sep = ":")
cluster_index <- match(dat$samples$lonlat, unique(dat$samples$lonlat))
first_index <- mapply(function(x) which(cluster_index == x)[1], 1:max(cluster_index))

# get mean PC value per cluster
pc_mean <- mapply(function(y) {
  mapply(mean, split(y, f = cluster_index))
}, split(pca$x[,1:4], f = col(pca$x[,1:4])))
colnames(pc_mean) <- c("a) PC1", "b) PC2", "c) PC3", "d) PC4")

# make plotting dataframe
plot_df <- subset(dat$samples, select = c("long", "lat"))[first_index,]
plot_df <- cbind(plot_df, pc_mean)

# get into long format for facetted plot
plot_df_long <- tidyr::gather(data = plot_df, key = component, value = value, 3:6, factor_key = TRUE)

# -----------------------------------
# PLOTTING

# define plotting parameters
col_country <- grey(0.95)
col_country_border <- grey(0.5)
size_country_border <- 0.5
col_sea <- grey(1.0)
resolution <- "coarse"
col_limits <- c(-4,4)
point_size <- 2.5
stroke_col <- grey(0.5)
stroke_size <- 0.25
col_vec <- rev(colorRampPalette(bobfunctions2::col_tim())(100))

# load country shapefiles
world_map <- getMap(resolution = resolution)

# basic map plot
plot_base <- ggplot() + theme_bw() + theme(panel.background = element_rect(fill = col_sea),
                                           panel.grid.major = element_blank(),
                                           panel.grid.minor = element_blank(),
                                           strip.background = element_blank(),
                                           strip.text = element_text(angle = 0, hjust = 0, size = 12))

# add country borders
plot_base <- plot_base + geom_polygon(aes(long, lat, group = group),
                                      size = size_country_border, color = col_country_border,
                                      fill = col_country, data = world_map)


# make separate inset plots
inset_list <- list()
for (i in 1:4) {
  
  # get colours for this principal component
  plot_df$col <- plot_df[,2+i]
  
  # create inset plot from base
  inset_list[[i]] <- plot_base + coord_cartesian(xlim = c(-5, 3), ylim = c(4,11)) +
    geom_point(aes(x = long, y = lat, fill = col), color = stroke_col, shape = 21,
               stroke = stroke_size, size = point_size, data = plot_df) +
    scale_fill_gradientn(colours = col_vec, limits = col_limits, guide = FALSE) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    theme(plot.margin = grid::unit(c(1,1,0,0), "mm"))
}

# create main facetted plot
plot1 <- plot_base + coord_cartesian(xlim = c(5, 33), ylim = c(-13,8)) +
  geom_point(aes(x = long, y = lat, fill = value), color = stroke_col, shape = 21,
             stroke = stroke_size, size = point_size, data = plot_df_long) +
  scale_fill_gradientn(colours = col_vec, name = "PC value", limits = col_limits) +
  xlab("longitude") + ylab("latitude") +
  facet_wrap(~component)

# add rect for inset
plot_combined <- plot1 + geom_rect(aes(xmin = x0, xmax = x1, ymin = y0, ymax = y1), color = "black",
                                   fill = NA, data = data.frame(x0 = 4, x1 = 16, y0 = -1, y1 = 8.6))

# add inset plots
for (i in 1:4) {
  plot_combined <- plot_combined + gg_inset(ggplotGrob(inset_list[[i]]),
                                            data = data.frame(component = names(plot_df)[2+i]),
                                            xmin = 4, xmax = 16, ymin = -1, ymax = 8.6)
}

# save to file
ggsave("figure3_PCA_maps/figure3_PCA_maps.pdf", plot = plot_combined, device = "pdf", width = 9, height = 7)
ggsave("figure3_PCA_maps/figure3_PCA_maps.png", plot = plot_combined, device = "png", width = 9, height = 7, dpi = 100)
