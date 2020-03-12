
# generate_figure1.R
#
# Author: Bob Verity
# Date: 2019-12-17
#
# Purpose:
# Produce 2D and 3D PCA plots. 3D version uses bespoke code for plotting 3D scatter using ggplot2.
#
# ------------------------------------------------------------------

# NOTE - uncomment these lines to install packages as needed
#devtools::install_github("mrc-ide/MIPAnalyzer", ref = "version1.0.0")
#devtools::install_github("bobverity/bobfunctions2", ref = "version1.0.0")

# load packages
library(devtools)
library(MIPanalyzer)
library(bobfunctions2)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)

# ------------------------------------------------------------------

# read in raw data
dat <- readRDS("source_data/biallelic_distances.rds")

# PCA on WSAFs
wsaf_impute <- get_wsaf(dat, impute = TRUE, FUN = mean)
pca <- pca_wsaf(wsaf_impute)

# plotting arguments
pointsize_2D <- 1
pointsize_3D <- 0.3
panel_letter_size <- 8

# create plotting dataframe
plot_df <- data.frame(x = pca$x[,1], y = pca$x[,2], z = pca$x[,3], country = dat$samples$Country)

# change order of countries. Changes order in which points are plotted so that DRC is behind
w <- match(plot_df$country, c("DRC", "Ghana", "Tanzania", "Uganda", "Zambia"))
plot_df <- plot_df[order(w),]

# define colour legend
col_pal <- c(grey(0.8), RColorBrewer::brewer.pal(8, "Dark2")[c(2,1,3,4)])

# make axis labels
x_lab <- sprintf("PC1 (%s%%)", round(pca$var[1], digits = 2))
y_lab <- sprintf("PC2 (%s%%)", round(pca$var[2], digits = 2))
z_lab <- sprintf("PC3\n(%s%%)", round(pca$var[3], digits = 2))

# produce 2D scatterplot
plot1 <- ggplot() + theme_bw(base_size = 8) +
  geom_point(aes(x = x, y = y, col = country), data = plot_df, size = pointsize_2D) + 
  xlab(x_lab) + ylab(y_lab) + ggtitle("a)") + theme(plot.title = element_text(size = panel_letter_size)) +
  scale_color_manual(values = col_pal, name = "Country", guide = FALSE)

# produce 3D scatterplot
plot2 <- gg3d_scatterplot(plot_df$x, plot_df$y, plot_df$z, colour = plot_df$country, size = pointsize_3D, theta = 135,
                          d = 2.0, x_lab = x_lab, y_lab = y_lab, z_lab = z_lab,
                          axis_lab_dist = 1.5, axis_lab_size = 1.8, grid_size = 0.15, zero_line_size = 0.4) +
  scale_color_manual(values = col_pal, name = "Country") +
  ggtitle("b)") + theme(plot.title = element_text(size = panel_letter_size),
                        legend.text = element_text(size = 8),
                        legend.title = element_text(size = 8))

# plot side-by-side
plot3 <- gridExtra::grid.arrange(plot1, plot2, ncol = 2, widths = c(2.1,3))

plot(plot3)

# save to file
file_ext <- c("eps", "pdf", "png")
for (i in seq_along(file_ext)) {
  ggsave(sprintf("figure1_PCA/figure1_PCA.%s", file_ext[i]),
         plot = plot3, device = file_ext[i], width = 179, height = 80, units = "mm")
}
