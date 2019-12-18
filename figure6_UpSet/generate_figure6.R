
# generate_figure6.R
#
# Author: Bob Verity
# Date: 2019-12-17
#
# Purpose:
# Generate plots of drug resistance combinations using the UpSet package
#
# ------------------------------------------------------------------

# NOTE - uncomment these lines to install packages as needed
#devtools::install_github("mrc-ide/MIPAnalyzer", ref = "version1.0.0")

# load packages
library(devtools)
library(MIPanalyzer)
library(RColorBrewer)
library(UpSetR)
library(rworldmap)

# ------------------------------------------------------------------
# PRE-PROCESS

# read in data
dat <- readRDS("source_data/dr_processed2.rds")

# subset to DRC only
dat <- filter_samples(dat, dat$samples$Country == "DRC")

# read in OJ monoclonals
coi_all <- readRDS("source_data/non_summariased_cois.rds")

# subset to chosen level K-means clustering, het thresholds etc.
coi_mono <- subset(coi_all, region_denom == 1 & gt == 0.1 & quantile0.975 == 1)

# subset MIP data to monoclonals
dat <- filter_samples(x = dat,
                      sample_filter = dat$samples$ID %in% coi_mono$name,
                      description = "subset to OJ monoclonals")

# add spatial jitter
n_samp <- nrow(dat$samples)
dat$samples$long_jitter <- dat$samples$long + rnorm(n_samp, sd = 0.1)
dat$samples$lat_jitter <- dat$samples$lat + rnorm(n_samp, sd = 0.1)

# ------------------------------------------------------------------
# DHPS UPSET

# get dhps mutants
get_mut_unique <- function(gene = "dhps") {
  mut_mat <- t(mapply(function(i) dat$samples$aa_seq[[i]][[gene]][,"s_mut"], 1:nrow(dat$samples)))
  mut_unique <- unique(as.vector(mut_mat))
  mut_unique[!is.na(mut_unique)]
}
mut_unique <- get_mut_unique("dhps")

# make upset matrix
get_upset_mat <- function(gene = "dhps", mut_unique) {
  ret <- t(mapply(function(i) mut_unique %in% dat$samples$aa_seq[[i]][[gene]][,"s_mut"], 1:nrow(dat$samples)))
  ret <- as.data.frame(ret*1)
  names(ret) <- mut_unique
  return(ret)
}
upset_mat <- get_upset_mat("dhps", mut_unique)

# account for G437A reference being resistant
upset_mat$G437A <- 1 - upset_mat$G437A
mut_unique[mut_unique == "G437A"] <- "A437G"
names(upset_mat) <- mut_unique

# account for situation of no mutations
upset_mat$none <- (rowSums(upset_mat) == 0)*1

# define plotting parameters
bp <- RColorBrewer::brewer.pal(11, "RdYlBu")
bar_cols <- c(bp[5], grey(0.5), bp[8], bp[4], bp[10], bp[2])

# save plot to file
pdf(file = "figure6_UpSet/figure6a_UpSet_dhps.pdf", onefile = FALSE, width = 4.5, height = 4)
UpSetR::upset(upset_mat, sets = c("none", mut_unique), keep.order = FALSE, main.bar.color = bar_cols, mb.ratio = c(0.5,0.5))
dev.off()

# ------------------------------------------------------------------
# DHPS MAP

# define map parameters
map_cols <- bar_cols[-2]
point_size <- 2.5
stroke_col <- grey(0.5)
stroke_size <- 0.25
col_sea <- "white"
resolution <- "coarse"
size_country_border <- 0.25
col_country_border <- grey(0.0)
col_country <- grey(0.9)

# define mutant sets for plotting map
mut_set <- NA
mut_set[upset_mat$A437G == 1 & rowSums(upset_mat) == 1] <- 1
mut_set[upset_mat$S436A == 1 & rowSums(upset_mat) == 1] <- 2
mut_set[upset_mat$A437G == 1 & upset_mat$K540E == 1 & rowSums(upset_mat) == 2] <- 3
mut_set[upset_mat$A437G == 1 & upset_mat$S436A == 1 & rowSums(upset_mat) == 2] <- 4
mut_set[rowSums(upset_mat) == 3] <- 5
dat$samples$mut_set <- mut_set

# drop NA rows
plot_df <- dat$samples[!is.na(dat$samples$mut_set),]

# load country shapefiles
world_map <- rworldmap::getMap(resolution = resolution)

# create basic map plot
plot_base <- ggplot() + theme_bw() +
  theme(panel.background = element_rect(fill = col_sea),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(angle = 0, hjust = 0, size = 12))

# create final map
plot1 <- plot_base +
  geom_polygon(aes(long, lat, group = group),
               size = size_country_border, color = col_country_border,
               fill = col_country, data = world_map) +
  geom_point(aes(x = long_jitter, y = lat_jitter, fill = as.factor(mut_set)), color = stroke_col,
             shape = 21, stroke = stroke_size, size = point_size, data = plot_df) +
  coord_cartesian(xlim = c(12, 33), ylim = c(-13,5)) +
  scale_fill_manual(values = map_cols, guide = FALSE) +
  xlab("longitude") + ylab("latitude") +
  ggtitle(expression(paste("b) ", italic("dhps"), " spatial distribution of mutant haplotypes")))

# save to file
ggsave("figure6_UpSet/figure6b_map_dhps.pdf", plot = plot1, device = "pdf", width = 4.5, height = 4)
ggsave("figure6_UpSet/figure6b_map_dhps.png", plot = plot1, device = "png", width = 4.5, height = 4, dpi = 100)

# ------------------------------------------------------------------
# CRT UPSET

# get crt mutants
mut_unique <- get_mut_unique("crt")

# make upset matrix
upset_mat <- get_upset_mat("crt", mut_unique)

# account for situation of no mutations
upset_mat$none <- (rowSums(upset_mat) == 0)*1

# define plotting parameters
bp <- RColorBrewer::brewer.pal(11, "RdYlBu")
bar_cols <- c(grey(0.5), "#CAB2D6", "#6A3D9A", bp[5], bp[4], bp[3], bp[8], bp[10], bp[2], bp[1])

# upset plot
pdf(file = "figure6_UpSet/figure6c_UpSet_crt.pdf", onefile = FALSE, width = 4.5, height = 4)
upset(upset_mat, sets = c("none", mut_unique), keep.order = FALSE, main.bar.color = bar_cols, mb.ratio = c(0.5,0.5))
dev.off()

# ------------------------------------------------------------------
# CRT MAP

# define map parameters
map_cols <- bar_cols[-1]

# define mutant sets for plotting map
mut_set <- NA
mut_set[upset_mat$A220S == 1 & rowSums(upset_mat) == 1] <- 1
mut_set[upset_mat$G32D == 1 & rowSums(upset_mat) == 1] <- 2
mut_set[upset_mat$M74I == 1 & rowSums(upset_mat) == 3] <- 3
mut_set[upset_mat$M74I == 1 & upset_mat$A220S & rowSums(upset_mat) == 4] <- 4
mut_set[upset_mat$M74I == 1 & upset_mat$N326S & rowSums(upset_mat) == 4] <- 5
mut_set[upset_mat$M74I == 1 & upset_mat$I356T & rowSums(upset_mat) == 4] <- 6
mut_set[upset_mat$M74I == 1 & upset_mat$A220S & upset_mat$I356T] <- 7
mut_set[upset_mat$M74I == 1 & upset_mat$A220S & upset_mat$N326S] <- 8
mut_set[upset_mat$M74I == 1 & upset_mat$A220S & upset_mat$N326S & upset_mat$F325C] <- 9
dat$samples$mut_set <- mut_set

# drop NA rows
plot_df <- dat$samples[!is.na(dat$samples$mut_set),]

# create final map
plot2 <- plot_base +
  geom_polygon(aes(long, lat, group = group),
               size = size_country_border, color = col_country_border,
               fill = col_country, data = world_map) +
  geom_point(aes(x = long_jitter, y = lat_jitter, fill = as.factor(mut_set)), color = stroke_col,
             shape = 21, stroke = stroke_size, size = point_size, data = plot_df) +
  coord_cartesian(xlim = c(12, 33), ylim = c(-13,5)) +
  scale_fill_manual(values = map_cols, guide = FALSE) +
  xlab("longitude") + ylab("latitude") +
  ggtitle(expression(paste("d) ", italic("crt"), " spatial distribution of mutant haplotypes")))
plot2

# save to file
ggsave("figure6_UpSet/figure6d_map_crt.pdf", plot = plot1, device = "pdf", width = 4.5, height = 4)
ggsave("figure6_UpSet/figure6d_map_crt.png", plot = plot1, device = "png", width = 4.5, height = 4, dpi = 100)
