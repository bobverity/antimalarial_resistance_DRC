
# generate_figureS7.R
#
# Author: Bob Verity
# Date: 2019-12-17
#
# Purpose:
# Plot within-sample allele frequencies of highly related samples.
#
# ------------------------------------------------------------------

library(ggplot2)

# ------------------------------------------------------------------

# read in data
dat <- readRDS("source_data/biallelic_distances.rds")

# plotting parameters
text_size <- 2.5

# subset samples to DRC only
w <- which(dat$samples$Country == "DRC")
dat$distance$spatial <- dat$distance$spatial[w,w]
dat$distance$inbreeding_dominant <- dat$distance$inbreeding_dominant[w,w]
dat <- filter_samples(dat, dat$samples$Country == "DRC")

# get within-sample allele frequencies
wsaf <- get_wsaf(dat)

# find highly related pairs
f_thresh <- 0.9
w <- which(dat$distance$inbreeding_dominant > f_thresh, arr.ind = TRUE)
