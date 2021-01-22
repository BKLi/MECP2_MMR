# -*- coding: utf-8 -*-
"""
Created on Fri Aug 21 10:46:37 2020

@author: bingkun
@project: MMR

* enrichment of epigenetic markers on weak/strong enhancers and negative elements
"""

#%%
import pandas as pd
import numpy as np
import seaborn as sns
from pathlib import Path
import re
from functools import reduce
from scipy import stats
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib as mpl
mpl.rcParams['figure.dpi']= 300

date="0118"

#%%
pre_process = False
#%%
if pre_process:
# read in & group element
        element_all = pd.read_csv(r'C:\Users\libin\UCSF\MMR\enrichment\ATAC_Elements.csv', sep=",")
        # group
        element_strong = element_all[element_all["enhancer"] == "strong"]
        element_weak = element_all[element_all["enhancer"] == "weak"]
        element_negative = element_all[element_all["enhancer"] == "negative"]
        
        element_strong = element_strong[["chr", "Start", "End"]]
        element_weak = element_weak[["chr", "Start", "End"]]
        element_negative = element_negative[["chr", "Start", "End"]]
        
        element_strong.to_csv(r'C:\Users\libin\UCSF\MMR\enrichment\element_strong', sep="\t", header=False, index=False)
        element_weak.to_csv(r'C:\Users\libin\UCSF\MMR\enrichment\element_weak', sep="\t", header=False, index=False)
        element_negative.to_csv(r'C:\Users\libin\UCSF\MMR\enrichment\element_negative', sep="\t", header=False, index=False)

#%%
input_dir = Path(r'C:\Users\libin\UCSF\MMR\enrichment')
file_paths = input_dir.glob('**/*pooled*.tab')
file_path_list = [str(path) for path in file_paths]
print (file_path_list)       
        
#%%
df_list = []
for file in file_path_list:
        file_name = file.split("\\")[-1]
        basename = "".join(re.findall(r'(.+)_summary', file_name))
        dataframe = pd.read_csv(file, sep="\t", skiprows=1, names=['chr','start','end','ChIP_H3K27ac_norm_pooled','Cut_Tag_CTCF_pooled','Cut_Tag_H3K27ac_pooled','Cut_Tag_H3K4me3_pooled','LM059_LM063'])
        #dataframe = dataframe.drop(['ChIP_H3K27ac_pooled', 'Cut_Tag_CTCF_pooled', 'Cut_Tag_H3K27ac_pooled', 'Cut_Tag_H3K4me3_pooled'], axis=1)
        #dataframe = dataframe.rename(columns={"KJ207": "CTCF_1", "KJ209": "CTCF_2", "KJ208": "H3K27ac_1", "KJ210": "H3K27ac_2", "KJ220": "H3K27ac_3", "KJ222": "H3K27ac_4", "LMA93": "ATAC_1", "LMA97": "ATAC_2", "NE009": "H3K4me3_1", "NE013": "H3K4me3_2"})
        #dataframe = dataframe[["CTCF_1", "CTCF_2", "H3K27ac_1", "H3K27ac_2", "H3K27ac_3", "H3K27ac_4", "H3K4me3_1", "H3K4me3_2", "ATAC_1", "ATAC_2"]]
        dataframe = dataframe.rename(columns={"LM059_LM063": "ATAC", "ChIP_H3K27ac_norm_pooled": "ChIP_H3K27ac", "Cut_Tag_H3K27ac_pooled": "Cut_Tag_H3K27ac", "Cut_Tag_CTCF_pooled": "CTCF", "Cut_Tag_H3K4me3_pooled": "H3K4me3"})
        dataframe = dataframe[['ATAC','ChIP_H3K27ac','Cut_Tag_H3K27ac','CTCF','H3K4me3']]
        print (dataframe.shape[0])
        dataframe_melt = pd.melt(dataframe, value_name="marker_signal", var_name="")
        dataframe_melt["type"] = basename

        ax = plt.figure(figsize=(16,10))
        plt.ylabel('signal', fontsize=15,fontname="Arial")
        plt.xlabel('{}'.format(basename), fontsize=15,fontname="Arial")

        plt.xticks(rotation=45, fontname="Arial", fontsize=12)
        plt.yticks(fontname="Arial", fontsize=12)
        ax = sns.lineplot(data=dataframe)
        #plt.savefig(r'C:\Users\libin\UCSF\MMR\enrichment/element_enrichment_0827_linePlot_{}.pdf'.format(basename), transparent=True) 
        
        df_list.append(dataframe_melt)
        dataframe_T = dataframe.T
        dataframe_T_melt = pd.melt(dataframe_T, value_name="marker_signal", var_name="")

        #dataframe["type"] = basename
        
#%%

dataframe_melt_concat = pd.concat(df_list)
dataframe_melt_concat.to_csv(r'C:\Users\libin\UCSF\MMR\enrichment\dataframe_melt_concat_{}.csv'.format(date), sep=",", index=False, header=False)

#%%
ax = plt.figure(figsize=(16,10))
plt.ylabel('', fontsize=12,fontname="Arial")
plt.xticks(rotation=45, fontname="Arial", fontsize=12)
plt.yticks(fontname="Arial", fontsize=12)
# sns.swarmplot(x="", y=y_ax, data=df_gr1_melt, color="blue", size=5, zorder=0)
# ax=sns.boxplot(x="", y=y_ax, data=df_gr1_melt, palette=["gray", "whitesmoke"], showfliers = False, width=[0.3], boxprops=dict(edgecolor='black',linewidth=3))
ax = sns.boxplot(x="", y="marker_signal", hue="type", data=dataframe_melt_concat, showfliers = False, palette='Pastel1')
plt.setp(ax.lines, color="grey", linewidth=0.5)
plt.setp(ax.spines.values(), color="black", linewidth=0.5)
ax = sns.swarmplot(x="", y="marker_signal", hue="type", data=dataframe_melt_concat, color='grey', size=1.5, dodge=True)

plt.savefig(r'C:\Users\libin\UCSF\MMR\enrichment/element_enrichment_{}.pdf'.format(date), transparent=True) 

#%%
## one plot for each marker
marker_list = ['ATAC','ChIP_H3K27ac','Cut_Tag_H3K27ac','CTCF','H3K4me3']
for m in marker_list:      
        dataframe_melt_concat_subset = dataframe_melt_concat[dataframe_melt_concat[""] == m]
        ax = plt.figure(figsize=(16,10))
        plt.ylabel('', fontsize=12,fontname="Arial")
        plt.xticks(rotation=45, fontname="Arial", fontsize=12)
        plt.yticks(fontname="Arial", fontsize=12)
        ax = sns.boxplot(x="", y="marker_signal", hue="type", data=dataframe_melt_concat_subset, showfliers = False, palette='Pastel1')
        plt.setp(ax.lines, color="grey", linewidth=0.5)
        plt.setp(ax.spines.values(), color="black", linewidth=0.5)
        ax = sns.swarmplot(x="", y="marker_signal", hue="type", data=dataframe_melt_concat_subset, color='grey', size=1.5, dodge=True)
        
        plt.savefig(r'C:\Users\libin\UCSF\MMR\enrichment/element_enrichment_{}_{}.pdf'.format(m, date), transparent=True) 

#%%
# Kruskal–Wallis test
stat_dic = {}
for m in marker_list:      
        dataframe_melt_concat_subset = dataframe_melt_concat[dataframe_melt_concat[""] == m]
        marker = "".join(dataframe_melt_concat_subset[""].to_list()[0])        
        dataframe_melt_concat_subset_control = dataframe_melt_concat_subset[dataframe_melt_concat_subset["type"] == "element_control"]
        dataframe_melt_concat_subset_neg = dataframe_melt_concat_subset[dataframe_melt_concat_subset["type"] == "element_negative"]
        dataframe_melt_concat_subset_strong = dataframe_melt_concat_subset[dataframe_melt_concat_subset["type"] == "element_strong"]
        dataframe_melt_concat_subset_weak = dataframe_melt_concat_subset[dataframe_melt_concat_subset["type"] == "element_weak"]
        
        marker_dic={}
        control_neg = stats.kruskal(dataframe_melt_concat_subset_control["marker_signal"], dataframe_melt_concat_subset_neg["marker_signal"], nan_policy="omit")[1]
        control_strong = stats.kruskal(dataframe_melt_concat_subset_control["marker_signal"], dataframe_melt_concat_subset_strong["marker_signal"], nan_policy="omit")[1]
        control_weak = stats.kruskal(dataframe_melt_concat_subset_control["marker_signal"], dataframe_melt_concat_subset_weak["marker_signal"], nan_policy="omit")[1]
        neg_strong = stats.kruskal(dataframe_melt_concat_subset_neg["marker_signal"], dataframe_melt_concat_subset_strong["marker_signal"], nan_policy="omit")[1]
        neg_weak = stats.kruskal(dataframe_melt_concat_subset_neg["marker_signal"], dataframe_melt_concat_subset_weak["marker_signal"], nan_policy="omit")[1]
        strong_weak = stats.kruskal(dataframe_melt_concat_subset_strong["marker_signal"], dataframe_melt_concat_subset_weak["marker_signal"], nan_policy="omit")[1]
        marker_dic["control_neg"] = control_neg
        marker_dic["control_strong"] = control_strong
        marker_dic["control_weak"] = control_weak
        marker_dic["neg_strong"] = neg_strong
        marker_dic["neg_weak"] = neg_weak
        marker_dic["strong_weak"] = strong_weak
        stat_dic[marker] = marker_dic
stat_dic_df = pd.DataFrame.from_dict(stat_dic,orient='index').reset_index().rename(columns={"index":"marker"})       


#%% summary of conservaton score
input_dir = Path(r'C:\Users\libin\UCSF\MMR\enrichment')
file_paths = input_dir.glob('*_cons.tab')
file_path_list = [str(path) for path in file_paths]
print (file_path_list)       
        
#%%
df_list = []
for file in file_path_list:
        file_name = file.split("\\")[-1]
        basename = "".join(re.findall(r'(.+)_summary_cons', file_name))
        dataframe = pd.read_csv(file, sep="\t", skiprows=1, names=['chr','start','end','cons_score'])
        #dataframe = dataframe.drop(['ChIP_H3K27ac_pooled', 'Cut_Tag_CTCF_pooled', 'Cut_Tag_H3K27ac_pooled', 'Cut_Tag_H3K4me3_pooled'], axis=1)
        #dataframe = dataframe.rename(columns={"KJ207": "CTCF_1", "KJ209": "CTCF_2", "KJ208": "H3K27ac_1", "KJ210": "H3K27ac_2", "KJ220": "H3K27ac_3", "KJ222": "H3K27ac_4", "LMA93": "ATAC_1", "LMA97": "ATAC_2", "NE009": "H3K4me3_1", "NE013": "H3K4me3_2"})
        #dataframe = dataframe[["CTCF_1", "CTCF_2", "H3K27ac_1", "H3K27ac_2", "H3K27ac_3", "H3K27ac_4", "H3K4me3_1", "H3K4me3_2", "ATAC_1", "ATAC_2"]]
        #dataframe = dataframe.rename(columns={"ATAC_pooled": "ATAC", "ChIP_H3K27ac_norm_pooled": "ChIP_H3K27ac", "Cut_Tag_H3K27ac_pooled": "Cut_Tag_H3K27ac", "Cut_Tag_CTCF_pooled": "CTCF", "Cut_Tag_H3K4me3_pooled": "H3K4me3"})
        dataframe = dataframe[['cons_score']]
        dataframe["type"] = basename
        #plt.savefig(r'C:\Users\libin\UCSF\MMR\enrichment/element_enrichment_0827_linePlot_{}.pdf'.format(basename), transparent=True) 
       
        df_list.append(dataframe)
dataframe_concat = pd.concat(df_list)
#%%        
ax = plt.figure(figsize=(12,12))
plt.ylabel('conservation_score', fontsize=15,fontname="Arial")
plt.xticks(rotation=45, fontname="Arial", fontsize=12)
plt.yticks(fontname="Arial", fontsize=12)
ax = sns.boxplot(x="type", y="cons_score", data=dataframe_concat, showfliers = False, palette='Pastel1')
ax = sns.swarmplot(x="type", y="cons_score", data=dataframe_concat, color='grey', size=1.5, dodge=True)

plt.savefig(r'C:\Users\libin\UCSF\MMR\enrichment/element_conservation_score.pdf', transparent=True) 

#%% summary of conservaton score -- max
input_dir = Path(r'C:\Users\libin\UCSF\MMR\enrichment')
file_paths = input_dir.glob('*max.out')
file_path_list = [str(path) for path in file_paths]
print (file_path_list)   

   
df_list = []
for file in file_path_list:
        file_name = file.split("\\")[-1]
        basename = "".join(re.findall(r'(.+).max.out', file_name))
        dataframe = pd.read_csv(file, delim_whitespace=True, names=['chr1','start1','end1','chr2', 'start2', 'end2', 'id', 'score'])
        dataframe = dataframe[['score']]
        dataframe["type"] = basename
        dataframe = dataframe[["type", "score"]]
        df_list.append(dataframe)
dataframe_concat = pd.concat(df_list)
dataframe_concat.to_csv(r'C:\Users\libin\UCSF\MMR\enrichment\conservation_concat.csv', sep=",", index=False, header=False)


ax = plt.figure(figsize=(12,12))
plt.ylabel('conservation_score', fontsize=15,fontname="Arial")
plt.xticks(rotation=45, fontname="Arial", fontsize=12)
plt.yticks(fontname="Arial", fontsize=12)
ax = sns.boxplot(x="type", y="score", data=dataframe_concat, showfliers = False, palette='Pastel1')
ax = sns.swarmplot(x="type", y="score", data=dataframe_concat, color='grey', size=1.5, dodge=True)

plt.savefig(r'C:\Users\libin\UCSF\MMR\enrichment/element_conservation_score_max.pdf', transparent=True) 
