
# fit_PrevMap.R
#
# Author: Bob Verity
# Date: 2019-12-17
#
# Purpose:
# Read in drug resistance data. Fit maps of allele frequencies over space using
# PrevMap. Save output to file.
#
# ------------------------------------------------------------------

# NOTE - uncomment these lines to install packages as needed
#devtools::install_github("mrc-ide/MIPAnalyzer", ref = "version1.0.0")

# load packages
library(devtools)
library(MIPanalyzer)
library(dplyr)
library(PrevMap)
library(splancs)
library(raster)
library(magrittr)

# ------------------------------------------------------------------

# set seed
set.seed(1)

# read in data
dat <- readRDS("source_data/dr_processed2.rds")

# subset to DRC
dat <- filter_samples(dat, dat$samples$Country == "DRC")

# extract mutation status in crt:K76T and dhps:K540E
mut_status <- mapply(function(i) {
  aa_df <- as.data.frame(do.call(rbind, dat$samples$aa_seq[[i]]))
  K76 <- as.character(subset(aa_df, gene == "crt" & codon_num_unique == 76)$s_mut)
  K540 <- as.character(subset(aa_df, gene == "dhps" & codon_num_unique == 540)$s_mut)
  return(c(K76, K540))
}, 1:nrow(dat$samples))

# get sample info into convenient dataframe
df_samp <- data.frame(lat = dat$samples$lat,
                      long = dat$samples$long,
                      clust = dat$samples$hv001,
                      K76 = mut_status[1,],
                      K540 = mut_status[2,])

# subset to samples with spatial information
df_samp <- subset(df_samp, !is.na(lat) & !is.na(long))

# group over clusters
df_clust <- df_samp %>%
  dplyr::group_by(clust) %>%
  dplyr::summarise(lat = lat[1],
                   long = long[1],
                   N = length(K76),
                   K76 = sum(!is.na(K76)),
                   K540 = sum(!is.na(K540))) %>%
  as.data.frame()

# logit-transform mutation prevalences
eps <- 0.5
df_clust$K76_logit <- log((df_clust$K76 + eps)/(df_clust$N - df_clust$K76 + eps))
df_clust$K540_logit <- log((df_clust$K540 + eps)/(df_clust$N - df_clust$K540 + eps))

# make list of model fits by maximum likelihood
fit <- list()
fit$K76T_fit <- PrevMap::linear.model.MLE(formula = K76_logit ~ 1, coords = ~long + lat,
                                          data = df_clust, start.cov.pars = c(3,3), kappa = 1)
fit$K540E_fit <- linear.model.MLE(formula = K540_logit ~ 1, coords = ~long + lat,
                                  data = df_clust, start.cov.pars = c(3,3), kappa = 1)

# create grid of points to map
poly <- cbind(c(12.2,12,32,32,18,12.4), c(-5.5,6,6,-14,-14,-6.2))
grid_pred <- splancs::gridpts(poly, xs = 0.25, ys = 0.25)
colnames(grid_pred) <- c("long","lat")

# make list of model predictions from fits
pred <- list()
for (i in 1:length(fit)) {
  
  # make model predictions
  pred_raw <- spatial.pred.linear.MLE(fit[[i]], grid_pred, scale.predictions = "prevalence",
                                      n.sim.prev = 1e3, standard.errors = TRUE)
  
  # make predictions into raster then dataframe
  rast <- raster::rasterFromXYZ(cbind(pred_raw$grid.pred, pred_raw$prevalence$predictions))
  df_rast <- as.data.frame(raster::rasterToPoints(rast))
  
  # make smooth raster (for use when producing contour lines)
  rast_smooth <- raster::disaggregate(raster::aggregate(rast, fact = 3), fact = 3, method = "bilinear")
  df_rast_smooth <- as.data.frame(raster::rasterToPoints(rast_smooth))
  
  # store in pred list
  pred[[i]] <- list(df_rast = df_rast,
                    df_rast_smooth = df_rast_smooth)
}

# save to file
saveRDS(pred, file = "generated_data/dr_PrevMap_predictions.rds")
