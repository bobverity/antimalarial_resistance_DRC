
# generate_figureS5.R
#
# Author: Bob Verity
# Date: 2019-12-17
#
# Purpose:
# Barplot of variance explained by each pricipal component.
#
# ------------------------------------------------------------------

# NOTE - uncomment these lines to install packages as needed
#devtools::install_github("mrc-ide/MIPAnalyzer", ref = "version1.0.0")

# load packages
library(devtools)
library(MIPanalyzer)

# ------------------------------------------------------------------

# read in data
dat <- readRDS("source_data/biallelic_distances.rds")

# PCA on WSAFs
wsaf_impute <- get_wsaf(dat, impute = TRUE, FUN = mean)
pca <- pca_wsaf(wsaf_impute)

# create plotting dataframe
n <- 20
df_plot <- data.frame(component = 1:n, variance = pca$var[1:n])

# produce plot
plot1 <- ggplot() + theme_bw() +
  geom_col(aes(x = component, y = variance), data = df_plot) +
  ylab("percent variance explained")

# save to file
ggsave("figureS5_PCA_variance/figureS5_PCA_variance.pdf", plot = plot1, device = "pdf", width = 5, height = 4)
ggsave("figureS5_PCA_variance/figureS5_PCA_variance.png", plot = plot1, device = "png", width = 5, height = 4, dpi = 100)