
# generate_figure4.R
#
# Author: Bob Verity
# Date: 2019-12-17
#
# Purpose:
# Generate histogram of estimated IBD.
#
# ------------------------------------------------------------------

# NOTE - uncomment these lines to install packages as needed
#devtools::install_github("mrc-ide/MIPAnalyzer", ref = "version1.0.0")
#devtools::install_github("bobverity/bobfunctions2", ref = "version1.0.0")

# load packages
library(devtools)
library(MIPanalyzer)
library(bobfunctions2)

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

# plotting parameters
fill_col <- grey(0.85)
inset_x <- c(0.5, 1)
inset_y <- c(20, 70)
breaks <- seq(0, 1, 0.02)
df_bracket <- data.frame(x0 = c(0.5,1,0.75,0.5), x1 = c(0.5,1,0.75,1),
                         y0 = c(2,2,5,5), y1 = c(5,5,20,5))

# create main histogram
plot1 <- ggplot() + theme_bw() +
  geom_histogram(aes(x = x, y = (..count..)/sum(..count..)*100),
                 breaks = breaks, fill = fill_col, color = "black",
                 data = data.frame(x = f_vec)) +
  xlab("IBD") + ylab("frequency (%)") +
  geom_segment(aes(x = x0, xend = x1, y = y0, yend = y1), data = df_bracket)

# create inset plot
plot2 <- ggplot() + theme_bw() +
  geom_histogram(aes(x = x, y = (..count..)/sum(..count..)*100),
                 breaks = breaks, fill = fill_col, color = "black",
                 data = data.frame(x = f_vec)) +
  coord_cartesian(xlim = c(0.5,1), ylim = c(0,0.002)) +
  xlab("IBD") + ylab("frequency (%)")

# create combined plot
plot_combined <- plot1 +
  bobfunctions2::gg_inset(ggplotGrob(plot2),
                          data = data.frame(z = 1),
                          xmin = inset_x[1], xmax = inset_x[2], ymin = inset_y[1], ymax = inset_y[2]) +
  geom_rect(aes(xmin = x0, xmax = x1, ymin = y0, ymax = y1), color = "black", fill = NA,
            data = data.frame(x0 = inset_x[1], x1 = inset_x[2], y0 = inset_y[1], y1 = inset_y[2]))

# save to file
ggsave("figure4_IBD_histogram/figure4_IBD_histogram.pdf", plot = plot_combined,
       device = "pdf", width = 6, height = 4.5)
ggsave("figure4_IBD_histogram/figure4_IBD_histogram.png", plot = plot_combined,
       device = "png", width = 6, height = 4.5, dpi = 100)
