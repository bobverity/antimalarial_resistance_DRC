
# generate_figure2.R
#
# Author: Bob Verity
# Date: 2019-12-17
#
# Purpose:
# Plot locus contributions for first few principal components of PCA.
#
# ------------------------------------------------------------------

# NOTE - uncomment these lines to install packages as needed
#devtools::install_github("mrc-ide/MIPAnalyzer", ref = "version1.0.0")

# load packages
library(devtools)
library(MIPanalyzer)
library(RColorBrewer)

# ------------------------------------------------------------------

# read in raw data
dat <- readRDS("source_data/biallelic_distances.rds")

# PCA on WSAFs
wsaf_impute <- get_wsaf(dat, impute = TRUE, FUN = mean)
pca <- pca_wsaf(wsaf_impute)

# get chromosome lengths
df_chrom <- Pf_chrom_lengths()
df_chrom$offset <- cumsum(df_chrom$length) - df_chrom$length
total_length <- sum(df_chrom$length)

# make plotting dataframe
comp_name_vec <- paste("Principal Component ", 1:4, sep = "")
df_plot <- data.frame(x = df_chrom$offset[dat$loci$CHROM_NUMERIC] + dat$loci$POS,
                      y = as.vector(pca$contribution[,1:4]),
                      chrom = df_chrom$chrom[dat$loci$CHROM_NUMERIC],
                      comp = rep(comp_name_vec, each = nrow(pca$contribution)),
                      geo = c("non-informative", "informative")[as.numeric(dat$loci$GEO)])

# choose point colours
point_cols <- c("#E31A1C", "#1F78B4")

# create basic plot
plot1 <- ggplot() +
  theme(panel.grid.minor.x = element_blank()) +
  geom_point(aes(x = x, y = y, color = geo), size = 0.7, data = df_plot) +
  facet_wrap(~comp, ncol = 1) +
  theme(strip.background = element_blank(), strip.text = element_text(size = 11, hjust = 0)) +
  xlab("SNP position (ordered chromosomes)") + ylab("relative contribution (%)") +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  theme(axis.ticks.x=element_blank(), axis.text.x=element_blank()) +
  scale_colour_manual(values = point_cols, name = "geographic\nstatus of site", guide = FALSE) +
  scale_x_continuous(limits = c(0,total_length), expand = c(0,0), breaks = cumsum(df_chrom$length)) +
  scale_y_continuous(limits = c(0,3), expand = c(0,0))

# create dataframe for annotating plot
df_ann <- data.frame(comp = comp_name_vec[c(2,4)],
                     lab = c("pfcrt", "dhps"),
                     x = c(7021158, 8667514),
                     y = 2.5)

# annotate plot
plot1 <- plot1 + geom_text(aes(x = x, y = y, label = lab), size = 3, fontface = 2, data = df_ann)

# save to file
ggsave("figure2_PCA_contribution/figure2_PCA_contribution.pdf", plot = plot1, device = "pdf",
       width = 12, height = 6)
ggsave("figure2_PCA_contribution/figure2_PCA_contribution.png", plot = plot1, device = "png",
       width = 12, height = 6, dpi = 100)
