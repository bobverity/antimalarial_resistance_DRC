
# generate_figureS1.R
#
# Author: Bob Verity
# Date: 2019-12-17
#
# Purpose:
# Plot depth (coverage) distributions for both MIP panels.
#
# ------------------------------------------------------------------

# NOTE - uncomment these lines to install packages as needed
#devtools::install_github("mrc-ide/MIPAnalyzer", ref = "version1.0.0")

# load packages
library(devtools)
library(MIPanalyzer)
library(ggplot2)
library(gridExtra)

# ------------------------------------------------------------------
# PRE-PROCESS

# read in data
dat1 <- readRDS("source_data/biallelic_processed0.rds")
dat2 <- readRDS("source_data/dr_processed0.rds")

# get coverage
c1 <- as.vector(dat1$coverage)
c2 <- as.vector(dat2$coverage)

# get quantiles
q1 <- quantile(c1, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
q2 <- quantile(c2, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)

# create summary vectors
summary1 <- c(sprintf("mean = %s", round(mean(c1, na.rm = TRUE)), digits = 2),
              sprintf("Q25 = %s", q1[1]),
              sprintf("Q50 (median) = %s", q1[2]),
              sprintf("Q75 = %s", q1[3]))
summary2 <- c(sprintf("mean = %s", round(mean(c2, na.rm = TRUE)), digits = 2),
              sprintf("Q25 = %s", q2[1]),
              sprintf("Q50 (median) = %s", q2[2]),
              sprintf("Q75 = %s", q2[3]))

# histogram of genome-wide panel coverage
plot1 <- ggplot() + theme_bw(base_size = 8) +
  geom_histogram(aes(x = log(x)/log(10)), fill = grey(0.5), color = grey(0.3),
                                bins = 20, data = data.frame(x = c1)) +
  xlim(c(0,4)) + ylim(c(0,3e5)) +
  xlab(expression(log[10]*"(depth)")) + ggtitle("genome-wide panel") +
  ggplot2::annotate("text", hjust = 1, x = 4, y = 3e5*c(1,0.9,0.8,0.7), label = summary1, size = 3)

# histogram of DR panel coverage
plot2 <- ggplot() + theme_bw(base_size = 8) +
  geom_histogram(aes(x = log(x)/log(10)), fill = grey(0.5), color = grey(0.3),
                 bins = 20, data = data.frame(x = c2)) +
  xlim(c(0,5)) + ylim(c(0,5e5)) +
  xlab(expression(log[10]*"(depth)")) + ggtitle("drug resistance panel") +
  ggplot2::annotate("text", hjust = 1, x = 5, y = 5e5*c(1,0.9,0.8,0.7), label = summary2, size = 3)

# create combined plot
plot_combined <- gridExtra::grid.arrange(plot1, plot2, ncol = 2)

# save to file
file_ext <- c("eps", "pdf", "png")
for (i in seq_along(file_ext)) {
  ggsave(sprintf("figureS1_coverage/figureS1_coverage.%s", file_ext[i]),
         plot = plot_combined, device = file_ext[i], width = 179, height = 80, units = "mm")
}
