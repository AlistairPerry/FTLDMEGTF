#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 13:25:03 2023

@author: alistairperry
"""

"""
Setup

"""


import os
import pandas as pd
import numpy as np
import scipy.io

import seaborn as sns
import matplotlib.pyplot as plt

import ptitprince as pt

# Single condition waveforms
LFPmat_data = scipy.io.loadmat("/Users/alistairperry/Documents/Cambridge/Project_2/Fix_T1cor/TF/B_Data/commconpatplac_wTGB_LFPcond_MMN_allROIs_wsem_con.mat")

# MMN waveforms
MMNmat_data = scipy.io.loadmat("/Users/alistairperry/Documents/Cambridge/Project_2/Fix_T1cor/TF/B_Data/commconpatplac_wTGB_LFPdiffcalcs_MMN_RIFG_wsem_wPatSubgrps.mat")

# Mean Df
df_MMN = pd.read_csv("/Users/alistairperry/Documents/Cambridge/Project_2/Fix_T1cor/TF/B_Data/commconpatplac_wTGB_LFPs_MMNdiffmean_ConsandPats_forJASP_fullsample.txt", delimiter="\t")

# Remove doubles
df_MMN = df_MMN.loc[df_MMN['is_double']==0] 

# Group colours
pal = ["#3498db", "#E69F00", "#00BA38"]

# OutDir
OutDir = "/Users/alistairperry/Documents/Cambridge/Project_2/Fix_T1cor/TF/C_Plots/"



"""
Single-condition waveforms - all ROIs

"""

ROI_labs = ["L IFG", "L STG", "L AUD", "R IFG", "R STD", "R AUD"]

# Create time vector
time_pts = np.arange(0, 251)


# Init plot
plt.figure(figsize=(20, 10))

for ind, roi in enumerate(ROI_labs):
    
    # Init subplot
    ax = plt.subplot(2, 3, ind + 1)
    
    # Extract and plot mean data
    m1_std = LFPmat_data["m1std_con_RIFG"][:,:,ind].reshape(251)
    m1_dev = LFPmat_data["m1dev_con_RIFG"][:,:,ind].reshape(251)

    ax.plot(time_pts, m1_std, color="blue", linewidth=1.5, label="rep6")
    ax.plot(time_pts, m1_dev, color="red", linewidth=1.5, label="dev")
    
    # Now SEM
    s1_std = LFPmat_data["s1std_con_RIFG"][ind,:].reshape(251)
    s1_dev = LFPmat_data["s1dev_con_RIFG"][ind,:].reshape(251)
    
    ax.fill_between(time_pts, m1_std-s1_std, m1_std+s1_std, color="blue", alpha=.15)
    ax.fill_between(time_pts, m1_dev-s1_dev, m1_dev+s1_dev, color="red", alpha=.15)
    
    # Stylise labels
    plt.title(ROI_labs[ind], fontsize=24, fontweight="bold")
    
    plt.xlabel("Time (ms)", fontsize=24, fontweight="bold")
    plt.ylabel("Amplitude (a.u)", fontsize=24, fontweight="bold")

    plt.yticks(fontsize=20, fontweight="bold")
    plt.xticks(ticks=[0, 50, 100, 150, 200, 250], labels=['-100', '0', '100', '200', '300', '400'], fontsize=20, fontweight="bold")

    plt.tight_layout()
    
    # Legend
    if ind == 0:
        plt.legend(fontsize=20)

# Save out
plt.savefig(os.path.join(OutDir, "LFP_bothconds_allROIS_onlycons.png"), dpi=300)
plt.show()


"""
MMN waveforms - per subgroup

"""

# Create time vector
time_pts = np.arange(0, 251)

# Init plot
fig, ax = plt.subplots(figsize=(12, 8))

grp_name = ["con", "pat_bv", "pat_psp"]
grp_lab = ["CON", "bvFTD", "PSP"]

for ind, grp in enumerate(grp_name):
    m1 = MMNmat_data[f"m1diff_{grp}_RIFG"]
    m1 = m1.reshape(251)
    
    s1 = MMNmat_data[f"s1diff_{grp}_RIFG"]
    s1 = s1.reshape(251)

    ax.plot(time_pts, m1, color=pal[ind], linewidth=3, label=grp_lab[ind])
    ax.fill_between(time_pts, m1-s1, m1+s1, color=pal[ind], alpha=.15)

# Legend
legend_properties = {'weight':'bold', 'size': 20}
plt.legend(loc=2, prop=legend_properties)

# Stylise labels
plt.xlabel("Time (s)", fontsize=24, fontweight="bold")
plt.ylabel("Mismatch response\n(rep6-dev)", fontsize=24, fontweight="bold")

plt.yticks(fontsize=20, fontweight="bold")
plt.xticks(ticks=[0, 50, 100, 150, 200, 250], labels=['-0.1', '0', '0.1', '0.2', '0.3', '0.4'], fontsize=20, fontweight="bold")

# Other
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)
plt.tight_layout()

# Save
plt.savefig(os.path.join(OutDir, "Mismatch_3grps_RIFG.tif"), dpi=300)
plt.show()


"""
Mean MMN

"""

# Init plot
f, ax = plt.subplots(figsize=(10, 8))


#Add box plots
ax=sns.boxplot(x="Diag", y="meanMMNcol_4", data=df_MMN, order=[3, 1, 2], width=.4, linewidth=3,
               zorder=10, showcaps=True, showfliers=False, whiskerprops={'linewidth':3, "zorder":10},
               saturation=1, orient="v")

# Set colours
ax = sns.set_palette(pal)

#Add individual data points with jitter
ax=sns.stripplot(x="Diag", y="meanMMNcol_4", data=df_MMN, order=[3, 1, 2], edgecolor="white",
size=15, jitter=1, zorder=0, orient="v")


# adding transparency to colors
for patch in ax.patches:
 r, g, b, a = patch.get_facecolor()
 patch.set_facecolor((r, g, b, .7))



#Manually change labels

plt.ylabel("Mean MMN (125-175ms)", fontsize=24, fontweight="bold")

plt.xlabel('Group', fontsize=24, fontweight="bold")

plt.xticks(ticks=[0, 1, 2], labels=["CON", "bvFTD", "PSP"], fontsize=20, fontweight="bold")

plt.yticks(fontsize=20, fontweight="bold")

ax.spines.right.set_visible(False)

ax.spines.top.set_visible(False)

plt.tight_layout()

plt.savefig(os.path.join(OutDir, "MMN_3grps_RIFG.tif"), dpi=300)

plt.show()