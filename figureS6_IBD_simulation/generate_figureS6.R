
# generate_figureS6.R
#
# Author: Bob Verity
# Date: 2019-12-17
#
# Purpose:
# Plot results of simulation-based analysis of identity by descent (IBD).
# Compares "true" (simulated) IBD against estimated value using maximum
# likelihood estimator.
#
# ------------------------------------------------------------------

library(ggplot2)

# ------------------------------------------------------------------

# load summary data from file
summary_df <- readRDS("source_data/sim_summary.rds")

# subset dataframe
plot_df <- subset(summary_df, type %in% c("truth", "estimated",
                                          "estimated_sub1", "estimated_sub3", "estimated_sub5"))

# rename some values
plot_df$m_text <- c("m = 0.00", "m = 0.25", "m = 0.50", "m = 1.00")[match(plot_df$m, c(0, 0.25, 0.5, 1))]
plot_df$mean_coi_text <- c("COI = 1", "COI = 2", "COI = 3")[match(plot_df$mean_coi, c(1e-12, 1.594, 2.822))]
plot_df$type_text <- c("'truth'", "MLE 1079SNPs", "MLE 500SNPs", "MLE 100SNPs", "MLE 20SNPs")[match(plot_df$type, c("truth", "estimated", "estimated_sub1", "estimated_sub3", "estimated_sub5"))]
plot_df$type_text <- factor(plot_df$type_text, levels = c("'truth'", "MLE 1079SNPs", "MLE 500SNPs", "MLE 100SNPs", "MLE 20SNPs"))

# produce plot of average IBD
plot1 <- ggplot(data = plot_df, aes(x = N, y = IBD_mean, color = type_text)) + theme_bw() +
  geom_line() +
  scale_x_continuous(trans ='log10') +
  ylim(c(0,1)) + 
  ylab("mean IBD") +
  scale_color_discrete(name = "IBD type") +
  facet_grid(vars(mean_coi_text), vars(m_text))

# save to file
ggsave("figureS6_IBD_simulation/figureS6_IBD_simulation.pdf", plot = plot1,
       device = "pdf", width = 8, height = 6)
ggsave("figureS6_IBD_simulation/figureS6_IBD_simulation.png", plot = plot1,
       device = "png", width = 8, height = 6, dpi = 100)
