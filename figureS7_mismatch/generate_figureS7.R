
# generate_figureS7.R
#
# Author: Bob Verity
# Date: 2019-12-17
#
# Purpose:
# Plot within-sample allele frequencies of highly related samples.
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

# loop through highly related pairs. Produce plot for each
plot_list <- list()
for (i in 1:nrow(w)) {
  
  # get plotting data
  df_plot <- data.frame(x = wsaf[w[i,1],], y = wsaf[w[i,2],])
  
  # count matches and mismatches
  counts <- rep(NA,4)
  counts[1] <- sum((df_plot$x < 0.5) & (df_plot$y >= 0.5))
  counts[2] <- sum((df_plot$x >= 0.5) & (df_plot$y >= 0.5))
  counts[3] <- sum((df_plot$x < 0.5) & (df_plot$y < 0.5))
  counts[4] <- sum((df_plot$x >= 0.5) & (df_plot$y < 0.5))
  pc <- round(sum(counts[2:3])/sum(counts)*100, digits = 1)
  
  # create basic plot
  plot1 <- ggplot() + theme_bw() +
    geom_point(aes(x = x, y = y), alpha = 0.2, data = df_plot) +
    xlab("REF-WSAF1") + ylab("REF-WSAF2") +
    geom_hline(aes(yintercept = 0.5), linetype = "dashed") +
    geom_vline(aes(xintercept = 0.5), linetype = "dashed")
  
  # annotate with counts
  plot1 <- plot1 + annotate("text", label = sprintf("mis-match\nn = %s", counts[1]),
                            size = text_size, color = "red", x = 0.25, y = 0.75) +
    annotate("text", label = sprintf("match\nn = %s", counts[2]),
             size = text_size, color = "blue", x = 0.75, y = 0.75) +
    annotate("text", label = sprintf("match\nn = %s", counts[3]),
             size = text_size, color = "blue", x = 0.25, y = 0.25) +
    annotate("text", label = sprintf("mis-match\nn = %s", counts[4]),
             size = text_size, color = "red", x = 0.75, y = 0.25) +
    ggtitle(sprintf("sample %s vs. %s", w[i,1], w[i,2]))
  
  plot_list[[i]] <- plot1
}


# create combined plot
plot_combined <- grid.arrange(grobs = plot_list)

# save to file
ggsave("figureS7_mismatch/figureS7_mismatch.pdf", plot = plot_combined, device = "pdf", width = 8, height = 10)
ggsave("figureS7_mismatch/figureS7_mismatch.png", plot = plot_combined, device = "png", width = 8, height = 10, dpi = 100)
