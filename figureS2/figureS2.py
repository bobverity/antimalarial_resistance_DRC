#!/usr/bin/env python
# coding: utf-8

# ## Supplemental Figures 2a/b

# In[40]:


# RUN
import sys
sys.path.append("/opt/src")
import mip_functions as mip
import json
import subprocess
import gzip
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import pandas as pd
import seaborn as sns
import numpy as np
import warnings
warnings.filterwarnings('ignore')


# Load control strain variant calls from pf3k crosses project

# In[2]:


linenum = 0
vcf_file = "/opt/work/pf3k.controls.vcf.gz"
with gzip.open(vcf_file, ) as infile:
    for line in infile:
        if line.startswith(b"##"):
            linenum += 1
        else:
            break
pf3k_calls = pd.read_table(vcf_file, skiprows=linenum).drop_duplicates()


# In[3]:


pf3k_calls.drop(["QUAL", "ID", "FILTER", "INFO", "FORMAT"], inplace=True, axis=1)
pf3k_calls.set_index(list(pf3k_calls.columns[:4]), inplace=True)


# Load variant calls from MIP experiments

# In[4]:


linenum = 0
vcf_file = "/opt/work/variants.vcf.gz"
with gzip.open(vcf_file, ) as infile:
    for line in infile:
        if line.startswith(b"##"):
            linenum += 1
        else:
            break
variant_calls = pd.read_table(vcf_file, skiprows=linenum).drop_duplicates()


# Filter to variants with single vcf records

# In[5]:


variant_calls = variant_calls.groupby(["#CHROM", "POS", "REF", "ALT"]).filter(lambda a: len(a) == 1)


# In[6]:


variant_calls.drop(["QUAL", "ID", "FILTER", "INFO", "FORMAT"], inplace=True, axis=1)
variant_calls.set_index(list(variant_calls.columns[:4]), inplace=True)


# Get control sample names and filter variants to control samples

# In[7]:


control_sample_ids = [c for c in variant_calls.columns
                     if "-".join(c.split("-")[:-2]) in ["D6", "D10"]]


# In[8]:


control_calls = variant_calls.loc[:, control_sample_ids]


# Load targeted alleles

# In[9]:


targets = pd.read_csv("/opt/project_resources/target_alleles.csv")


# In[10]:


target_alleles = targets.set_index(["#CHROM", "POS", "REF", "ALT"])


# Filter pf3k calls to targeted alleles

# In[11]:


pf3k_calls = pf3k_calls.reindex(target_alleles.index, axis=0)


# Filter experimental calls to targeted alleles

# In[12]:


controls = control_calls.reindex(target_alleles.index).fillna(".")


# Define function to within sample allele frequency for MIP variants

# In[13]:


def get_freq(gt):
    vals = gt.split(":")
    gen = vals[0]
    if gen == ".":
        return np.nan
    else:
        gen = gen.split("/")
        ad = list(map(float, vals[1].split(",")))
        counts = dict(zip(gen, ad))
        # call variants only when coverage depth is >= 10
        if sum(ad) >= 10:
            try:
                return round(counts["1"]/sum(ad), 3)
            except KeyError:
                return 0
        else:
            return np.nan
    


# In[14]:


control_freq = controls.applymap(get_freq)


# Define strain mix ratios as determined by spectrophotometer

# In[15]:


d_mix = {"HB3": 0.14, "7G8": 0.13, "DD2": 0.06}


# exctract genotypes from pf3k calls

# In[16]:


strain_genotypes = pf3k_calls.dropna(how="any").applymap(lambda a: 1 if a.split(":")[0] == "1/1"
                                   else 0 if a.split(":")[0] == "0/0"
                                   else np.nan)


# Filter pf3k calls where any strain is missing the genotype call

# In[17]:


strain_genotypes = strain_genotypes.loc[~(strain_genotypes["DD2"].isnull() |
                    strain_genotypes["HB3"].isnull()|
                    strain_genotypes["7G8"].isnull())]


# Calculate expected frequency of each targeted SNP based on pf3k genotypes

# In[18]:


strain_genotypes["Expected"] = sum(d_mix[c] * strain_genotypes[c] for c in strain_genotypes.columns)


# Get the mean frequency of each SNP within the control mixtures (observed)

# In[19]:


d_average = control_freq.dropna(how="all", axis=0).T.apply(["mean", "std"], axis=0).T


# In[20]:


d_average.index.names = ["#CHROM", "POS", "REF", "ALT"]


# Merge expected and observed data

# In[21]:


d_average = d_average.merge(strain_genotypes[["Expected"]], left_on=["#CHROM", "POS", "REF", "ALT"],
              right_on=["#CHROM", "POS", "REF", "ALT"])


# Print correlation value between the expected and observed for SNPs with expected value > 0

# In[22]:


d_average.rename(columns={"mean": "Observed"}, inplace=True)


# In[23]:


print("Pearson's correlation value betweent observed and expected allele frequencies are ")
print(d_average.loc[d_average["Expected"] > 0, ["Observed", "Expected"]].corr())


# Save files

# In[25]:


control_freq.to_csv("/opt/analysis/controls_observed_WSAF.csv")
d_average.to_csv("/opt/analysis/controls_observed_average_and_expected_WSAF.csv")


# Plot and save Supplementary Figure 2a

# In[26]:


fig, ax = plt.subplots()
dat = d_average.reset_index().sort_values("Expected")
dat.index = range(len(dat))
sns.lineplot(data=dat[["Observed", "Expected"]],
            ax=ax)
fig.set_dpi(300)
ax.set_ylabel("Within Sample Allele Frequency")
ax.set_xlabel("SNPs")


# In[27]:


fig.savefig("/opt/analysis/figureS2a_observed_vs_expected_WSAF.png", frameon=True, dpi=300)


# Save plot data

# In[28]:


dat.to_csv("/opt/analysis/figureS2a_data.csv")


# Calculate relative error

# In[34]:


import copy


# In[37]:


d_filt = copy.deepcopy(d_average.loc[d_average["Expected"] > 0])


# In[38]:


d_filt["Relative Error"] = abs((d_filt["Observed"] - d_filt["Expected"]) / d_filt["Expected"])


# Plot relative error

# In[41]:


fig, ax = plt.subplots()
sns.barplot(x="Expected", y="Relative Error", data=d_filt, ci=95, ax=ax)
fig.set_dpi(300)
ax.set_xlabel("Expected Frequency")


# In[42]:


fig.savefig("/opt/analysis/figureS2b_relative_error_vs_expected_WSAF.png", frameon=True, dpi=300)


# Save data

# In[43]:


d_filt.to_csv("/opt/analysis/figureS2b_data.csv")

