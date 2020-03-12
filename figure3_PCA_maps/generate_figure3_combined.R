
# generate_figure3_combined.R
#
# Author: Bob Verity
# Date: 2019-12-17
#
# Purpose:
# Produce maps of first few principal components.
#
# ------------------------------------------------------------------

# load packages
library(gridExtra)

# load plot grobs from file
plot1 <- readRDS("figure3_PCA_maps/figure3_part1.rds")
plot2.1 <- readRDS("figure3_PCA_maps/figure3_part2_plot1.rds")
plot2.2 <- readRDS("figure3_PCA_maps/figure3_part2_plot2.rds")

# create combined plot
lay <- rbind(c(rep(1,7), rep(2,4)),
             c(rep(1,7), rep(3,4)))
plot_combined <- gridExtra::grid.arrange(plot1, plot2.1, plot2.2, layout_matrix = lay)

# save to file
file_ext <- c("eps", "pdf", "png")
for (i in seq_along(file_ext)) {
  ggsave(sprintf("figure3_PCA_maps/figure3_combined.%s", file_ext[i]),
         plot = plot_combined, device = file_ext[i], width = 179, height = 100, units = "mm")
}
