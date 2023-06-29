# get_Liza_dr_data.R
#
# Author: Bob Verity
# Date: 2023-06-29
#
# Inputs: Andria_DRdata/dr_data_v2.rds
#
# Outputs: Liza_DRdata/dr_data_Liza_v1.rds
#
# Purpose:
# Reads in pre-filtered drug resistance data, originally prepared for Andria.
# Wrangles into format for use by Liza.
#
# ------------------------------------------------------------------

library(dplyr)

# read in raw data
dat_raw <- readRDS("Andria_DRdata/dr_data_v2.rds")

# filter to K540E mutations
dat <- dat_raw |>
  select(lat, long, Year, codon_540) |>
  filter(!is.na(codon_540)) |>
  group_by(Year, lat, long) |>
  summarise(n_mut = sum(codon_540 == "K540E"),
            n_total = n()) |>
  ungroup()

# save to file
if (FALSE) {
  saveRDS(dat, "Liza_DRdata/dr_data_Liza_v1.rds")
}