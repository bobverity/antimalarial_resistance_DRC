# get_dr_data_v2.R
#
# Author: Bob Verity
# Date: 2023-04-04
#
# Purpose:
# Variation of script for extracting dhps drug resistance mutation data, in
# which rather than subsetting to monoclonals based on OJ's Real McCOIL analysis
# we just look for samples that are monoclonal at this particular gene. So
# keeping anything whether the dominant allele is at >90% relative frequency.
# This leads to an expanded sample set of >500 samples, compared with the ~100
# for RMC monoclonals.
#
# ------------------------------------------------------------------

# load packages
#devtools::install_github("mrc-ide/MIPanalyzer", ref = "v1.0.0")
library(MIPanalyzer)
#devtools::install_github("bobverity/bobfunctions2", ref = "4bb9b78a9e5ed3b084ca928785adbfa0caf3885a")
library(bobfunctions2)
library(tidyverse)

# read in raw data
dat <- readRDS("source_data/dr_processed.rds")

# filter loci to dhps locus
dat2 <- MIPanalyzer::filter_loci(dat, dat$loci$CHROM == "chr8" &
                                   dat$loci$POS >= 548200 &
                                   dat$loci$POS <= 550616)

# for each sample-locus combination, find out if counts are distributed such
# that any count makes up less than 90% of the total coverage
n_samp <- nrow(dat2$samples)
n_loc <- nrow(dat2$loci)
hint_polyclonal <- matrix(NA, n_samp, n_loc)
for (i in 1:n_samp) {
  for (j in 1:n_loc) {
    x <- dat2$counts[,i,j]
    p <- x / sum(x, na.rm = TRUE)
    hint_polyclonal[i,j] <- any(p < 0.9 & p > 0.1, na.rm = TRUE)
  }
}

# drop these loci
w <- which(hint_polyclonal, arr.ind = TRUE)
for (i in seq_len(nrow(w))) {
  dat2$counts[, w[i,1], w[i,2]] <- rep(NA, 4)
  dat2$coverage[w[i,1], w[i,2]] <- NA
}

# explore the effects of filtering loci and samples by coverage based on
# thresholds
MIPanalyzer::explore_filter_coverage_loci(dat2)
MIPanalyzer::explore_filter_coverage_samples(dat2)

# filter samples based on coverage
dat3 <- MIPanalyzer::filter_coverage_samples(dat2)

# get matrix of nucleotides for each sample-locus combination, and which are
# non-REF
n_samp <- nrow(dat3$samples)
n_loc <- nrow(dat3$loci)
nuc_mat <- matrix(NA, n_samp, n_loc)
nonREF_mat <- matrix(NA, n_samp, n_loc)
for (j in 1:n_loc) {
  # get vector of all alternative alleles at this locus
  loc_ref <- as.character(dat3$loci$REF[j])
  loc_alts <- strsplit(as.character(dat3$loci$ALT[j]), ",")[[1]]
  for (i in 1:n_samp) {
    x <- dat3$counts[,i,j]
    if (!all(is.na(x))) {
      w <- which.max(x)
      if (w == 1) {
        nuc_mat[i,j] <- loc_ref
        nonREF_mat[i,j] <- FALSE
      } else {
        nuc_mat[i,j] <- loc_alts[w - 1]
        nonREF_mat[i,j] <- TRUE
      }
    }
  }
}

# filter loci to only those with any mutants in this sample
any_mutants <- apply(nonREF_mat, 2, function(x) any(x, na.rm = TRUE))
nuc_mat <- nuc_mat[,any_mutants]
nonREF_mat <- nonREF_mat[,any_mutants]
dat4 <- MIPanalyzer::filter_loci(dat3, any_mutants)

# read in list of information for genes of interest. including dhps
gene_list <- readRDS("source_data/gene_list.rds")

# get the relative position of our final loci from the start of the dhps gene
# (first position = 1)
rel_pos <- dat4$loci$POS - gene_list$dhps$start + 1

# define which codons we're interested in (we can find this out from the codons
# at rel_pos, but hard-coded here for simplicity)
codon_focus <- c(436, 437, 540, 581, 613)

# find the DNA positions corresponding to these codons
pos_focus <- matrix(which(gene_list$dhps$codon_num %in% codon_focus),
                    ncol = 3, byrow = TRUE)

# get the amino acids of the REF sequence at each codon
aa_ref <- mapply(function(i) {
  codon_seq <- substr(gene_list$dhps$dna_seq, pos_focus[i,1], pos_focus[i,3])
  bobfunctions2::dna_to_aa(codon_seq, 2)
}, 1:nrow(pos_focus))

# store the mutation name for all samples
mut_name_mat <- matrix(NA, nrow(nuc_mat), length(codon_focus))
for (i in 1:nrow(nuc_mat)) {
  
  # make a new complete gene sequence with mutations
  dna_this <- gene_list$dhps$dna_seq
  codon_mask <- rep(FALSE, length(codon_focus))
  for (j in 1:ncol(nuc_mat)) {
    if (is.na(nuc_mat[i,j])) {
      codon_mask[match(gene_list$dhps$codon_num[rel_pos[j]], codon_focus)] <- TRUE
    } else {
      substr(dna_this, rel_pos[j], rel_pos[j]) <- nuc_mat[i,j]
    }
  }
  
  # get the amino acid sequence at each codon
  aa_this <- mapply(function(j) {
    codon_seq <- substr(dna_this, pos_focus[j,1], pos_focus[j,3])
    bobfunctions2::dna_to_aa(codon_seq, 2)
  }, 1:nrow(pos_focus))
  
  # mask out NA nucleotide values
  aa_this[codon_mask] <- NA
  
  # name mutation by comparing against ref
  aa_compare <- (aa_this == aa_ref)
  mut_name <- rep(NA, length(codon_focus))
  if (any(aa_compare == TRUE, na.rm = TRUE)) {
    w <- which(aa_compare == TRUE)
    mut_name[w] <- "REF"
  }
  if (any(aa_compare == FALSE, na.rm = TRUE)) {
    w <- which(aa_compare == FALSE)
    mut_name[w] <- paste0(aa_ref[w], codon_focus[w], aa_this[w])
  }
  
  # store name
  mut_name_mat[i,] <- mut_name
}

# convert to named data.frame
mut_name_df <- as.data.frame(mut_name_mat)
names(mut_name_df) <- sprintf("codon_%s", codon_focus)

# append mutation name to sample information
dat_final <- dat4$samples %>%
  dplyr::select(ID, Country, ADM1NAME, DHSREGNA, Year, lat, long) %>%
  bind_cols(mut_name_df) %>%
  dplyr::filter(Country == "DRC")

# save to file
if (FALSE) {
  saveRDS(dat_final, file = "Andria_DRdata/dr_data_v2.rds")
}
