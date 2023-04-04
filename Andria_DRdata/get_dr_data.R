# get_dr_data.R
#
# Author: Bob Verity
# Date: 2023-04-03
#
# Purpose:
# Read in processed data and subset to monoclonals in the same way as in paper
# (using OJ's monoclonal calls). Make long form data.frame of sample information
# alongside drug resistance mutation information
#
# ------------------------------------------------------------------

# load packages
#devtools::install_github("mrc-ide/MIPanalyzer", ref = "v1.0.0")
library(MIPanalyzer)
library(tidyverse)

#### Read in data and ubset to monoclonals

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

# filter samples to columns we care about
df_samples <- dat$samples %>%
  dplyr::select(ID, Country, ADM1NAME, DHSREGNA, Year, lat, long)


#### Summarise DR mutations per sample

ret_list <- list()
for (i in 1:nrow(df_samples)) {
  df_mut <- do.call(rbind, dat$samples$aa_seq[[i]]) %>%
    as.data.frame()
  ret_list[[i]] <- expand_grid(df_samples[i,], df_mut)
}

# make final data.frame
df_ret <- ret_list %>%
  bind_rows()

# save to file
if (FALSE) {
  saveRDS(df_ret, file = "Andria_DRdata/dr_data.rds")
}