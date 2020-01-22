
# generate_figure5.R
#
# Author: Bob Verity
# Date: 2019-12-17
#
# Purpose:
# Generate plot of IBD vs spatial distance, and map of highly related links.
#
# ------------------------------------------------------------------

# NOTE - uncomment these lines to install packages as needed
#devtools::install_github("mrc-ide/MIPAnalyzer", ref = "version1.0.0")
#devtools::install_github("bobverity/bobfunctions2", ref = "version1.0.0")

# load packages
library(devtools)
library(MIPanalyzer)
library(dplyr)
library(magrittr)
library(bobfunctions2)
library(rgdal)
library(ggplot2)
library(rworldmap)

# ------------------------------------------------------------------

# read in data
dat <- readRDS("source_data/biallelic_distances.rds")

# subset to DRC
w <- which(dat$samples$Country == "DRC")
dat <- filter_samples(dat, 1:nrow(dat$samples) %in% w)
dat$distance$spatial <- dat$distance$spatial[w,w]
dat$distance$inbreeding_dominant <- dat$distance$inbreeding_dominant[w,w]

# get estimated relatedness as vector
f_vec <- as.vector(as.dist(t(dat$distance$inbreeding_dominant)))

# get average IBD within and between clusters
c <- dat$samples$hv001
u <- unique(c)
IBD_within <- IBD_between <- list()
for (i in 1:length(u)) {
  w1 <- which(c == u[i])
  z <- dat$distance$inbreeding_dominant[w1,w1]
  IBD_within[[i]] <- z[!is.na(z)]
  w2 <- which(c != u[i])
  z <- dat$distance$inbreeding_dominant[w1,w2]
  IBD_between[[i]] <- z[!is.na(z)]
}

# get mean and standard error of IBD split by spatial separation
d_vec <- as.vector(as.dist(t(dat$distance$spatial)))
breaks <- c(0, 1, seq(0,5000,100)[-1])
df_summary <- data.frame(f = f_vec,
                         d = d_vec,
                         c = cut(d_vec, breaks, include.lowest = TRUE)) %>%
  dplyr::group_by(c) %>%
  dplyr::summarise(f_mean = mean(f), f_SE = sd(f)/sqrt(length(f)))

# get midpoints of spatial intervals, and upper and lower 95% CIs on IBD
df_summary$mid <- midpoints(breaks)[match(df_summary$c, levels(df_summary$c))]
df_summary$f_UL <- df_summary$f_mean + 1.96*df_summary$f_SE
df_summary$f_LL <- df_summary$f_mean - 1.96*df_summary$f_SE

# create first plot
plot1 <- ggplot(data = df_summary) + theme_bw() +
  geom_point(aes(x = mid, y = f_mean)) +
  geom_segment(aes(x = mid, xend = mid, y = f_LL, yend = f_UL)) +
  xlab("spatial distance (km)") + ylab("mean IBD") +
  ggtitle("a)")

# ------------------------------------------------------------------
# Map of highly related links

# plotting parameters
col_country <- grey(0.95)
col_country_border <- grey(0.5)
size_country_border <- 0.5
col_sea <- grey(1.0)
shape_resolution <- "coarse"

# subset pairs to highly related
w <- which(dat$distance$inbreeding_dominant > 0.9, arr.ind = TRUE)
df_highrelated <- data.frame(x0 = dat$samples$long[w[,1]],
                             x1 = dat$samples$long[w[,2]],
                             y0 = dat$samples$lat[w[,1]],
                             y1 = dat$samples$lat[w[,2]],
                             f = dat$distance$inbreeding_dominant[w])

# get great circle distances and midpoints
df_highrelated$gc_dist <- bobfunctions2::lonlat_to_bearing(df_highrelated$x0,
                                                           df_highrelated$y0,
                                                           df_highrelated$x1,
                                                           df_highrelated$y1)$gc_dist
df_highrelated$x_mid <- (df_highrelated$x0 + df_highrelated$x1)/2
df_highrelated$y_mid <- (df_highrelated$y0 + df_highrelated$y1)/2

# get unique cluster locations, and vector of clusters containing clonal pairs
unique_clust <- unique(dat$samples$hv001)
w_clust <- cbind(dat$samples$hv001[w[,1]], dat$samples$hv001[w[,2]])
clonal_clust <- w_clust[w_clust[,1] == w_clust[,2],1]

# create dataframe of cluster locations. col = 1 if contains clonal pairs
df_clust <- dat$samples[match(unique_clust, dat$samples$hv001), c("long", "lat")]
df_clust$col <- 0
df_clust$col[match(clonal_clust, unique_clust)] <- 1
df_clust <- df_clust[order(df_clust$col),]

# get spatial distance between highly related pairs
d_hr <- dat$distance$spatial[w]

# load country shapefiles
world_map <- getMap(resolution = shape_resolution)

# read in shapefile of DRC water bodies
shape_raw <- rgdal::readOGR("source_data/DRCWater")
shape_raw@data$id <- rownames(shape_raw@data)
df_water <- merge(ggplot2::fortify(shape_raw, region = "id"), shape_raw@data, by = "id")

# create basic map plot
plot_base <- ggplot() + theme_bw() + theme(panel.background = element_rect(fill = col_sea),
                                           panel.grid.major = element_blank(),
                                           panel.grid.minor = element_blank(),
                                           strip.background = element_blank(),
                                           strip.text = element_text(angle = 0, hjust = 0, size = 12))

# create final map
plot2 <- plot_base +
  geom_polygon(aes(long, lat, group = group),
               size = size_country_border, color = col_country_border,
               fill = col_country, data = world_map) +
  coord_cartesian(xlim = c(12, 32), ylim = c(-13,5)) +
  geom_polygon(aes(x = long, y = lat, group = group), fill = grey(0), data = df_water) +
  geom_point(aes(x = long, y = lat), color = grey(0.7), data = df_clust) +
  geom_segment(aes(x = x0, xend = x1, y = y0, yend = y1, color = f), size = 1, data = df_highrelated) +
  xlab("longitude") + ylab("latitude") + ggtitle("b)") +
  scale_color_continuous(name = "relatedness (F)", limits = c(0.9, 1)) +
  geom_text(aes(x = x_mid, y = y_mid, label = round(gc_dist),
                hjust = c("left", "right", "left", "right", "left"),
                vjust = c("top", "top", "top", "bottom", "bottom")), size = 3.5,
            data = subset(df_highrelated, gc_dist > 1)) +
  annotate("text", x = 15, y = 3, label = "Congo river") +
  geom_segment(aes(x = 15, xend = 18.5, y = 2.5, yend = 1), arrow = arrow(length = unit(0.2, "cm")))


# create combined plot
plot_combined <- grid.arrange(plot1, plot2, nrow = 1, widths = c(1.3,2))

# save to file
ggsave("figure5_IBD_distance/figure5_IBD_distance.pdf", plot = plot_combined,
       device = "pdf", width = 10, height = 4.5)
ggsave("figure5_IBD_distance/figure5_IBD_distance.png", plot = plot_combined,
       device = "png", width = 10, height = 4.5, dpi = 100)
