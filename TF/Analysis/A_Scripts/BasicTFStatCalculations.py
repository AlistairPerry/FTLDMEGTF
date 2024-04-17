#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  7 13:05:40 2022

@author: alistairperry
"""


'''

Setup

'''

# Import libraries

import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import statannot

import ptitprince as pt



#Output Dir

OutDir = "/Users/alistairperry/Documents/Cambridge/Project_2/Fix_T1cor/TF/C_Plots/"


#Load Data

df = pd.read_csv("/Users/alistairperry/Documents/Cambridge/Project_2/Fix_T1cor/TF/B_Data/adv_ssst_newmaxf_fixICA_wfids_250_noTGBcon_FullDiffPowTable.csv")

df_bothreps = pd.read_csv("/Users/alistairperry/Documents/Cambridge/Project_2/Fix_T1cor/TF/B_Data/adv_ssst_newmaxf_fixICA_wfids_250_noTGBcon_FullPowCondTable_ConsPats.csv")

df_allreps = pd.read_csv("/Users/alistairperry/Documents/Cambridge/Project_2/Fix_T1cor/TF/B_Data/adv_ssst_newmaxf_fixICA_wfids_250_FullPowCondTable_ConsPats_allreps.csv")

df_beta2 = pd.read_csv("/Users/alistairperry/Documents/Cambridge/Project_2/Fix_T1cor/TF/B_Data/adv_ssst_newmaxf_fixICA_wfids_250_noTGBcon_FullDiffPowTable_hbeta.csv")

df_evo = pd.read_csv("/Users/alistairperry/Documents/Cambridge/Project_2/Fix_T1cor/TF/B_Data/adv_ssst_newmaxf_fixICA_wfids_250_EVOnew_noTGBcon_FullDiffPowTable.csv")


#Remove doubles straight away
df = df.loc[df['double']==0] 

df_allreps = df_allreps.loc[df_allreps['double']==0]


#Separate table for cons

df_cons = df.loc[df['Con_Pat']==1]

df_cons = df_cons.drop('ID', axis=1)

df_cons['ID'] = df_cons.index


df_bothreps_cons = df_bothreps.loc[df_bothreps['Con_Pat']==1]




regions = ["LIFG", "LSTG", "LAUD", "RIFG", "RSTG", "RAUD"]

reg_labs = ["L IFG", "L STG", "L AUD", "R IFG", "R STG", "R AUD"]


# Group colours

pal = ["#3498db", "#E69F00", "#00BA38"]


'''

Dep T Test [RIFG]

'''


#Robustness across control populations

from scipy.stats import ttest_ind

# tgb_con = df_cons.loc[df_cons['Study']==1]
# mem_con = df_cons.loc[df_cons['Study']==2]

# ttest_ind(tgb_con['RIFG'], mem_con['RIFG'])

# ax = sns.boxplot(x="Study", y="RIFG", data=df_cons, color="blue")

# ax = sns.swarmplot(x="Study", y="RIFG", data=df_cons, color=".25")

# # adding transparency to colors
# for patch in ax.patches:
#  r, g, b, a = patch.get_facecolor()
#  patch.set_facecolor((r, g, b, .7))
 
# plt.show()


#Paired T-Test


def PairedCondTest(df, regions):
    
    
    from scipy.stats import ttest_rel    

    
    t_all = []
    p_all = []
    
    
    for ind, reg in enumerate(regions):
        
        IV_std = ''.join([reg, "_Std"])
        IV_dev = ''.join([reg, "_Dev"])
        
        
        res = ttest_rel(df[IV_std], df[IV_dev])
        
        t_all.append(res[0])
        p_all.append(res[1])
        
    
    #Combine into output df    
    
    stats_df = np.stack((np.array(t_all), np.array(p_all)), axis=1)
    
    stats_output = pd.DataFrame(stats_df, index = regions, columns=["t", "p"])

    return stats_output
        
        
stats_output = PairedCondTest(df_bothreps_cons, regions)

stats_output.to_csv(path_or_buf=OutDir + "TF_PairedCondTests.csv", sep=',')

#RM ANOVA

from statsmodels.stats.anova import AnovaRM

df_cons_long = pd.melt(df_cons, id_vars='ID', value_vars=['LIFG', 'LSTG', 'LAUD', 'RIFG', 'RSTG', 'RAUD'])

df_cons_long = df_cons_long.rename(columns={"variable": "Region", "value": "TF"})


rm_model = AnovaRM(data=df_cons_long, depvar='TF',subject='ID', within=['Region']).fit()

print(rm_model.summary())


#Plot

ax = sns.boxplot(x="Region", y="TF", data=df_cons_long, color="blue")

ax = sns.swarmplot(x="Region", y="TF", data=df_cons_long, color=".25")


# adding transparency to colors
for patch in ax.patches:
 r, g, b, a = patch.get_facecolor()
 patch.set_facecolor((r, g, b, .7))


plt.axhline(y=0, color='r', linestyle='--')

plt.ylabel(r'$\Delta$ Beta Pow (Std-Dev) ', fontsize=12, fontweight="bold")

plt.show()


# Plot separately 

# plots are ok for now - just need to change 

def withincond_plot(region):
    
    plt.figure(figsize=(7,6))
    
    ax = sns.boxplot(y=region, data=df_cons, color="blue")

    ax = sns.swarmplot(y=region, data=df_cons, color=".25", size=15)

    # adding transparency to colors
    for patch in ax.patches:
     r, g, b, a = patch.get_facecolor()
     patch.set_facecolor((r, g, b, .3))


    plt.axhline(y=0, color='r', linestyle='--')
    
    plt.title(reg_labs[ind], fontsize=28, fontweight="bold")

    plt.ylabel(r'$\Delta$ Beta1 Power (rep6-dev)', fontsize=24, fontweight="bold")
    
    plt.yticks(weight = 'bold', size = 24)
    
    
    plt.savefig(os.path.join(OutDir, f"TFcon_conddiff_{region}.tif"))
    

for ind, reg in enumerate(regions):
    
    withincond_plot(reg)
    
    
# RIFG - box plot only
f, ax = plt.subplots(figsize=(10, 8))


#ax=pt.half_violinplot(x="RIFG", data=df_cons, palette=pal, bw=.5, cut=0., scale="area", width=.6, inner=None, orient="h")


ax=sns.boxplot(y="RIFG", data=df_cons, color="black", width=.15,linewidth=3, 
               zorder=10, showcaps=True, boxprops={'facecolor':'none', 
                                                   "zorder":10}, 
               showfliers=False, whiskerprops={'linewidth':3, "zorder":10},
               saturation=1, orient="v")

    
#Add individual data points with jitterÍ

ax=sns.stripplot(y="RIFG", data=df_cons, palette=pal, edgecolor="white",
size=17.5, jitter=1, zorder=0, orient="v")

# plt.ylim(-0.45, 0.45)

#Style labels
plt.title("R IFG", fontsize=32, fontweight="bold")

plt.ylabel(r'$\Delta$ Beta1 Power (rep6-dev) ', fontsize=26, fontweight="bold")

plt.yticks(fontsize=26, fontweight="bold")

plt.axhline(y=0, color='r', linestyle='--')

ax.spines.right.set_visible(False)

ax.spines.top.set_visible(False)

plt.tight_layout()

plt.savefig(os.path.join(OutDir, "TFcon_conddiff_RIFG_sing.tif"))
    

# All

find_max = np.max(abs(df_cons[regions]), axis=0).max()

# Init plot
plt.figure(figsize=(12.5, 7.5))

ax=sns.boxplot(x="Region", y="TF", data=df_cons_long, color="black", width=.4,linewidth=3, 
                   zorder=10, showcaps=True, boxprops={'facecolor':'none', 
                                                       "zorder":10}, 
                   showfliers=False, whiskerprops={'linewidth':3, "zorder":10},
                   saturation=1, orient="v")

        
#Add individual data points with jitterÍ
ax=sns.stripplot(x="Region", y="TF", data=df_cons_long, color=pal[0], edgecolor="white", size=15, jitter=1, zorder=0, orient="v")

#Style labels
plt.xlabel("Region", fontsize=24, fontweight="bold")
plt.ylabel(r'$\Delta$ Beta1 Power (rep6-dev) ', fontsize=24, fontweight="bold")
plt.yticks(fontsize=20, fontweight="bold")
plt.xticks(fontsize=20, fontweight="bold")
plt.axhline(y=0, color='r', linestyle='--')
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)
    
plt.tight_layout()

plt.savefig(os.path.join(OutDir, "TFcon_conddiff_beta1_all.png"), dpi=300)

'''

 Group comparisons

'''

import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import MultiComparison

import scikit_posthocs as sp


df_nodbl = df.loc[df['double']==0]


#ANOVA

# Fit and summarize ols model

model = ols('RIFG ~ C(Diag)',data=df_nodbl).fit()


# type-1 anova summary 

table_type_1_TF = sm.stats.anova_lm(model, typ=1)

table_type_1_TF.to_csv(path_or_buf=OutDir + "ANOVA_TF_3grps.csv", sep=',')


#Post-hoc tests

comparison_TF = MultiComparison(df_nodbl["RIFG"], df_nodbl["Diag"])
comparison_results_TF = comparison_TF.tukeyhsd()
print(comparison_results_TF.summary())


#Sidak

stats = sp.posthoc_ttest(df_nodbl, val_col="RIFG", group_col="Diag", p_adjust='sidak')

stats.to_csv(path_or_buf=OutDir + "ANOVA_TF_3grps_posthoctests.csv", sep=',')


#T-values (Con vs PSP)

stats_unc = ttest_ind(df_nodbl.loc[df_nodbl['Diag'] == 2, 'RIFG'], df_nodbl.loc[df_nodbl['Diag'] == 3, 'RIFG'])


#Print
f = open(OutDir + "ANOVA_TF_3grps_posthoctests_unc.txt", 'w')
print(f't= {stats_unc[0]}', file=f)
print(f'p= {stats_unc[1]}', file=f)
f.close()


#ANCOVA(S)

model_ancova = ols('RIFG ~ C(Diag) + Age',data=df_nodbl).fit()


# type-1 an(c)ova summary 

table_type_1_TF_ancova_Age = sm.stats.anova_lm(model_ancova, typ=1)

table_type_1_TF_ancova_Age.to_csv(path_or_buf=OutDir + "ANOVA_TF_3grps_AgeCov.csv", sep=',')


df_nodbl['RIFG_Dev_bl_log10'] = np.log10(df_nodbl['RIFG_Dev_bl'])

df_nodbl['RIFG_Std_bl_log10'] = np.log10(df_nodbl['RIFG_Std_bl'])


# Write out for BF stats
df_nodbl.to_csv(path_or_buf="/Users/alistairperry/Documents/Cambridge/Project_2/Fix_T1cor/TF/B_Data/adv_ssst_newmaxf_fixICA_wfids_250_noTGBcon_FullPowCondTable_ConsPats_plusBLlogpow.csv", index = False)


model_ancova = ols('RIFG ~ C(Diag) + RIFG_Dev_bl_log10', data=df_nodbl).fit()

table_type_1_TF_ancova_Devbl = sm.stats.anova_lm(model_ancova, typ=1)

table_type_1_TF_ancova_Devbl.to_csv(path_or_buf=OutDir + "ANOVA_TF_3grps_DevBLCov.csv", sep=',')



# Fit and summarize ols model

model = ols('RIFG_Dev_bl_log10 ~ C(Diag)',data=df_nodbl).fit()


# type-1 anova summary 

table_type_1_Dev_bl = sm.stats.anova_lm(model, typ=1)

table_type_1_Dev_bl.to_csv(path_or_buf=OutDir + "ANOVA_TF_Dev_bl_3grps.csv", sep=',')


stats = sp.posthoc_ttest(df_nodbl, val_col="RIFG_Dev_bl_log10", group_col="Diag", p_adjust='sidak')

stats.to_csv(path_or_buf=OutDir + "ANOVA_TF_Dev_bl_3grps_posthoctests.csv", sep=',')



# Plot

my_colors = ["#3498db", "#E69F00", "#00BA38"]

f, ax = plt.subplots(figsize=(10, 8))
ax = sns.boxplot(x="Diag", y="RIFG_Dev_bl_log10", data=df_nodbl, order=[3, 1, 2])

ax = sns.set_palette(my_colors)

ax = sns.swarmplot(x="Diag", y="RIFG_Dev_bl_log10", data=df_nodbl, order=[3, 1, 2], color=".25")


# adding transparency to colors
for patch in ax.patches:
 r, g, b, a = patch.get_facecolor()
 patch.set_facecolor((r, g, b, .7))
 
plt.xticks(ticks=[0, 1, 2], labels=["CON", "bvFTD", "PSP"])

plt.ylabel("Baseline dev Beta1 Power (log10)")
  
plt.show()



#Rainclouds


#Convert data to lists

conloc = df_nodbl.loc[df_nodbl['Diag']==3]
bvloc = df_nodbl.loc[df_nodbl['Diag']==1]
psploc = df_nodbl.loc[df_nodbl['Diag']==2]

a = list(conloc['RIFG'])
b = list(bvloc['RIFG'])
c = list(psploc['RIFG'])


#Generate indexes for two columns

size_a = len(a)
size_b = len(b)
size_c = len(c)

a_ind = list(np.tile(1, [size_a]))
b_ind = list(np.tile(2, [size_b]))
c_ind = list(np.tile(3, [size_c]))


#Transform data into a "long" structure
#Combine column data and indexes into a new pandas df

d = {'group': a_ind + b_ind + c_ind, 'score': a + b + c}

df_new = pd.DataFrame(d)

    
#Plot cloud distributions (half violin)
dx="group"; dy="score"

f, ax = plt.subplots(figsize=(9, 6))


ax=pt.half_violinplot(x=dx, y=dy, data=df_new, palette=pal, bw=.5, cut=0.,
                      scale="area", width=.6, inner=None, orient="v")


#Add box plots
ax=sns.boxplot(x=dx, y=dy, data=df_new, color="black", width=.15, linewidth=3,
               zorder=10, showcaps=True, boxprops={'facecolor':'none', 
                                                   "zorder":10}, 
               showfliers=False, whiskerprops={'linewidth':3, "zorder":10},
               saturation=1, orient="v")

    
#Add individual data points with jitter
ax=sns.stripplot(x=dx, y=dy, data=df_new, palette=pal, edgecolor="white",
size=10, jitter=1, zorder=0, orient="v")

# Add significance whiskers
statannot.add_stat_annotation(
    ax,
    data=df_new,
    x="group",
    y="score",
    box_pairs=[(3, 1), (3, 2), (1, 2)],
    test="t-test_ind",
    text_format="star",
    loc="inside",
)


#Manually change labels
plt.ylabel(r'$\Delta$ Beta1 Power (rep6-dev) ', fontsize=18, fontweight="bold")

plt.xlabel('Group', fontsize=18, fontweight="bold")

plt.xticks(ticks=[0, 1, 2], labels=["CON", "bvFTD", "PSP"], fontsize=18, fontweight="bold")

plt.yticks(fontsize=16, fontweight="bold")

ax.spines.right.set_visible(False)

ax.spines.top.set_visible(False)

plt.tight_layout()

plt.savefig(os.path.join(OutDir, "TFcon_3grps_RIFG.tif"), dpi=300)

plt.show()



"""
Beta Repetition Dynamics

"""

# Group to get condition means per group
df_grouped = (
    df_allreps[['Diagnosis', 'Condition','TF_scond']].groupby(['Diagnosis', 'Condition'])
    .agg(['mean'])
)

df_grouped = df_grouped.droplevel(axis=1, level=0).reset_index()

# Init plot
f, ax = plt.subplots(figsize=(12, 8))

# Line plot
sns.lineplot(data=df_grouped, x="Condition", y="mean", hue="Diagnosis", palette=pal, hue_order=[3, 1, 2], linewidth=3, marker="o", markersize=12)


# Manually change labels
plt.ylabel("Mean Beta1 Power (% change)", fontsize=24, fontweight="bold")
plt.xlabel('Repetition', fontsize=24, fontweight="bold")
plt.xticks(ticks=[1, 2, 3, 4, 5, 6, 7], labels=["Dev", "rep1", "rep2",  "rep3", "rep4", "rep5", "rep6"], fontsize=20, fontweight="bold")
plt.yticks(fontsize=20, fontweight="bold")

# Legend
legend_properties = {'weight':'bold', 'size': 20}
plt.legend(labels=["CON", "bvFTD", "PSP"], loc="best", prop=legend_properties)

# Other
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)

plt.tight_layout()

# Save out
plt.savefig(os.path.join(OutDir, "TF_allreps_3grps_py.tif"), dpi=300)
plt.show()



"""
Baseline Beta1 Power
"""

# Init plot
f, ax = plt.subplots(figsize=(10, 8))


#Add box plots
ax=sns.boxplot(x="Diag", y="RIFG_Dev_bl_log10", data=df_nodbl, order=[3, 1, 2], width=.4, linewidth=3,
               zorder=10, showcaps=True, showfliers=False, whiskerprops={'linewidth':3, "zorder":10},
               saturation=1, orient="v")

# Set colours
ax = sns.set_palette(pal)

#Add individual data points with jitter
ax=sns.stripplot(x="Diag", y="RIFG_Dev_bl_log10", data=df_nodbl, order=[3, 1, 2], edgecolor="white",
size=15, jitter=1, zorder=0, orient="v")


# adding transparency to colors
for patch in ax.patches:
 r, g, b, a = patch.get_facecolor()
 patch.set_facecolor((r, g, b, .7))

# Add significance whiskers
statannot.add_stat_annotation(
    ax,
    data=df_nodbl,
    x="Diag",
    y="RIFG_Dev_bl_log10",
    box_pairs=[(3, 1), (3, 2), (1, 2)],
    order=[3, 1, 2],
    test="t-test_ind",
    text_format="star",
    loc="inside",
)


# Manually change labels
plt.ylabel("Baseline Beta1 Power (log10) (Dev)", fontsize=24, fontweight="bold")
plt.xlabel('Group', fontsize=24, fontweight="bold")
plt.xticks(ticks=[0, 1, 2], labels=["CON", "bvFTD", "PSP"], fontsize=20, fontweight="bold")
plt.yticks(fontsize=20, fontweight="bold")

# Other
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)
plt.tight_layout()

# Save out
plt.savefig(os.path.join(OutDir, "TFcon_3grps_RIFG_bldev_log10.tif"), dpi=300)
plt.show()


"""
Beta2
"""

df_cons_beta2 = df_beta2.loc[df_beta2['Con_Pat']==1]

df_cons_beta2_long = pd.melt(df_cons_beta2, id_vars='ID', value_vars=['LIFG', 'LSTG', 'LAUD', 'RIFG', 'RSTG', 'RAUD'])

df_cons_beta2_long = df_cons_beta2_long.rename(columns={"variable": "Region", "value": "TF"})


# Init plot
plt.figure(figsize=(12.5, 7.5))

ax=sns.boxplot(x="Region", y="TF", data=df_cons_beta2_long, color="black", width=.4,linewidth=3, 
                   zorder=10, showcaps=True, boxprops={'facecolor':'none', 
                                                       "zorder":10}, 
                   showfliers=False, whiskerprops={'linewidth':3, "zorder":10},
                   saturation=1, orient="v")

        
#Add individual data points with jitterÍ
ax=sns.stripplot(x="Region", y="TF", data=df_cons_beta2_long, color=pal[0], edgecolor="white", size=15, jitter=1, zorder=0, orient="v")

#Style labels
plt.xlabel("Region", fontsize=24, fontweight="bold")
plt.ylabel(r'$\Delta$ Beta2 (22-30Hz) Power (rep6-dev)', fontsize=24, fontweight="bold")
plt.yticks(fontsize=20, fontweight="bold")
plt.xticks(fontsize=20, fontweight="bold")
plt.axhline(y=0, color='r', linestyle='--')
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)
    
plt.tight_layout()

plt.savefig(os.path.join(OutDir, "TFcon_conddiff_beta2_all.png"), dpi=300)


"""
Evoked
"""

df_cons_evo = df_evo.loc[df_evo['Con_Pat']==1]

# Single RIFG Plot
f, ax = plt.subplots(figsize=(10, 8))


#ax=pt.half_violinplot(x="RIFG", data=df_cons, palette=pal, bw=.5, cut=0., scale="area", width=.6, inner=None, orient="h")

ax=sns.boxplot(y="RIFG", data=df_cons_evo, color="black", width=.15,linewidth=3, 
               zorder=10, showcaps=True, boxprops={'facecolor':'none', 
                                                   "zorder":10}, 
               showfliers=False, whiskerprops={'linewidth':3, "zorder":10},
               saturation=1, orient="v")

    
#Add individual data points with jitterÍ

ax=sns.stripplot(y="RIFG", data=df_cons_evo, palette=pal, edgecolor="white",
size=17.5, jitter=1, zorder=0, orient="v")

# plt.ylim(-0.45, 0.45)

#Style labels
# plt.title("R IFG", fontsize=32, fontweight="bold")

plt.ylabel(r'$\Delta$ Evoked Beta1 Power (rep6-dev) ', fontsize=28, fontweight="bold")

plt.yticks(fontsize=28, fontweight="bold")

plt.axhline(y=0, color='r', linestyle='--')

ax.spines.right.set_visible(False)

ax.spines.top.set_visible(False)

plt.tight_layout()

plt.savefig(os.path.join(OutDir, "TFcon_conddiff_RIFG_sing_EVO.tif"))

# df_cons_evo_long = pd.melt(df_cons_evo, id_vars='ID', value_vars=['LIFG', 'LSTG', 'LAUD', 'RIFG', 'RSTG', 'RAUD'])

# df_cons_evo_long = df_cons_evo_long.rename(columns={"variable": "Region", "value": "TF"})


# # Init plot
# plt.figure(figsize=(12.5, 7.5))

# ax=sns.boxplot(x="Region", y="TF", data=df_cons_evo_long, color="black", width=.4,linewidth=3, 
#                    zorder=10, showcaps=True, boxprops={'facecolor':'none', 
#                                                        "zorder":10}, 
#                    showfliers=False, whiskerprops={'linewidth':3, "zorder":10},
#                    saturation=1, orient="v")

        
# #Add individual data points with jitterÍ
# ax=sns.stripplot(x="Region", y="TF", data=df_cons_evo_long, color=pal[0], edgecolor="white", size=15, jitter=1, zorder=0, orient="v")

# #Style labels
# plt.xlabel("Region", fontsize=24, fontweight="bold")
# plt.ylabel(r'$\Delta$ Evoked Beta1 Power (rep6-dev)', fontsize=24, fontweight="bold")
# plt.yticks(fontsize=20, fontweight="bold")
# plt.xticks(fontsize=20, fontweight="bold")
# plt.axhline(y=0, color='r', linestyle='--')
# ax.spines.right.set_visible(False)
# ax.spines.top.set_visible(False)
    
# plt.tight_layout()

