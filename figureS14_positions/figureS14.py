#!/usr/bin/env python
# coding: utf-8

# ## Supplemental Figures 2a/b

# In[1]:


# RUN
import sys
sys.path.append("/opt/src")
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import matplotlib.patches as patches
import pandas as pd
import seaborn as sns
import numpy as np
import warnings
warnings.filterwarnings('ignore')


# Load targeted alleles

# In[2]:


targets = pd.read_csv("/opt/project_resources/target_alleles.csv")


# Load chromosome sizes

# In[4]:


refs = pd.read_csv("/opt/work/chromosome_sizes.csv", index_col=0)


# In[5]:


mip_targets = targets.merge(refs)


# Label target types

# In[6]:


mip_targets["Type"] = np.nan
mip_targets.loc[mip_targets["Geo"], "Type"] = "Geographically Informative"
mip_targets["Type"].fillna("Neutral", inplace=True)


# Specify colors for target types

# In[7]:


mip_targets.loc[mip_targets["Geo"], "Color"] = "r"
mip_targets["Color"].fillna("b", inplace=True)


# Plot figure

# In[8]:


fig, ax = plt.subplots(7,2, sharex=True)
for i in range(1, 15):
    chrom = "chr" + str(i)
    coord = (((i - 1) // 2), ((i - 1) % 2), )
    data=mip_targets.loc[mip_targets["#CHROM"]==chrom]
    height = 0.5
    width = 1000
    ax[coord].bar(x="POS", data=data.loc[data["Geo"]],
                                              height=height, snap=False, width=width, color="Color",
                                              label="Geographically Informative")
    ax[coord].bar(x="POS", data=data.loc[~data["Geo"]],
                                              height=0.5, snap=False, width=width, color="Color",
                                              label="Neutral")
    ax[coord].set_facecolor('white')
    ax[coord].grid(False)
    if coord[0] != 6:
        ax[coord].tick_params(axis="both", bottom=False, top=False, left=False, right=False)
    else:
        ax[coord].tick_params(axis="y", bottom=False, top=False, left=False, right=False)
        ax[coord].tick_params(axis="x", which="minor", size=5)
        ax[coord].set_xlabel("Genomic Position")
    ax[coord].yaxis.set_ticklabels([])
    rect = patches.Rectangle((0, 0), data["chrom_size"].max(), height, zorder=0, 
                            edgecolor="none", facecolor="lightgrey")
    ax[coord].yaxis.set_ticklabels([])
    
    ax[coord].text(data["chrom_size"].max(), height/2, chrom)
    ax[coord].add_patch(rect)
fig.set_dpi(600)
fig.tight_layout()


# Save figure

# In[9]:


fig.savefig("/opt/analysis/figureS14_target_positions.png", frameon=True, dpi=300)


# Save data

# In[10]:


mip_targets.to_csv("/opt/analysis/figureS14_data.csv")

