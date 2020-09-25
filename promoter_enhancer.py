# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 20:15:27 2020

@ author: bingkun
@ summary plots for promoters acting as enhancers
"""

import pandas as pd
import seaborn as sns
from scipy import stats
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib as mpl
mpl.rcParams['figure.dpi']= 300

#%%

element_ann = pd.read_csv(r'C:\Users\libin\UCSF\MMR\manual_annotation_complete_update_refseq_0801.csv')
element_ann = element_ann[['chr','Start','End','annotation']]
#element_pmt_all = element_ann[element_ann["annotation"] == "promoter"]

#%%
element_strong_signal = pd.read_csv(r'C:\Users\libin\UCSF\MMR\enrichment\element_strong_summary_pooled.tab', delim_whitespace=True,
                                    skiprows=1, names=['chr','Start','End','ATAC_pooled','ChIP_H3K27ac_norm_pooled','Cut_Tag_CTCF_pooled','Cut_Tag_H3K27ac_pooled','Cut_Tag_H3K4me3_pooled'])
element_weak_signal = pd.read_csv(r'C:\Users\libin\UCSF\MMR\enrichment\element_weak_summary_pooled.tab', delim_whitespace=True,
                                    skiprows=1, names=['chr','Start','End','ATAC_pooled','ChIP_H3K27ac_norm_pooled','Cut_Tag_CTCF_pooled','Cut_Tag_H3K27ac_pooled','Cut_Tag_H3K4me3_pooled'])
element_neg_signal = pd.read_csv(r'C:\Users\libin\UCSF\MMR\enrichment\element_negative_summary_pooled.tab', delim_whitespace=True,
                                    skiprows=1, names=['chr','Start','End','ATAC_pooled','ChIP_H3K27ac_norm_pooled','Cut_Tag_CTCF_pooled','Cut_Tag_H3K27ac_pooled','Cut_Tag_H3K4me3_pooled'])
#%%
element_pos_signal = pd.concat([element_weak_signal, element_strong_signal], axis=0, ignore_index=True)
element_pos_ann = pd.merge(element_ann, element_pos_signal, on=['chr','Start','End'], how="inner")
element_pos_pmt = element_pos_ann[element_pos_ann["annotation"] == "promoter"]
element_pos_pmt = element_pos_pmt.rename(columns={"ATAC_pooled": "ATAC", "ChIP_H3K27ac_norm_pooled": "ChIP_H3K27ac", "Cut_Tag_H3K27ac_pooled": "Cut_Tag_H3K27ac", "Cut_Tag_CTCF_pooled": "CTCF", "Cut_Tag_H3K4me3_pooled": "H3K4me3"})
element_pos_pmt = element_pos_pmt[['ATAC','ChIP_H3K27ac','Cut_Tag_H3K27ac','CTCF','H3K4me3']]
element_pos_pmt_melt = pd.melt(element_pos_pmt, value_name="marker_signal", var_name="")
element_pos_pmt_melt["type"] = "positive"

element_neg_ann = pd.merge(element_ann, element_neg_signal, on=['chr','Start','End'], how="inner")
element_neg_pmt = element_neg_ann[element_neg_ann['annotation'] == 'promoter']
element_neg_pmt = element_neg_pmt.rename(columns={"ATAC_pooled": "ATAC", "ChIP_H3K27ac_norm_pooled": "ChIP_H3K27ac", "Cut_Tag_H3K27ac_pooled": "Cut_Tag_H3K27ac", "Cut_Tag_CTCF_pooled": "CTCF", "Cut_Tag_H3K4me3_pooled": "H3K4me3"})
element_neg_pmt = element_neg_pmt[['ATAC','ChIP_H3K27ac','Cut_Tag_H3K27ac','CTCF','H3K4me3']]
element_neg_pmt_sample = element_neg_pmt.sample(n=element_pos_pmt.shape[0])
element_neg_pmt_sample_melt = pd.melt(element_neg_pmt_sample, value_name="marker_signal", var_name="")
element_neg_pmt_sample_melt["type"] = "negative"

element_pmt_melt_all = pd.concat([element_pos_pmt_melt, element_neg_pmt_sample_melt])
element_pmt_melt_all.to_csv(r'C:\Users\libin\UCSF\MMR\enrichment\element_pmt_all.csv', sep=",", index=False, header=False)
#%%

x = plt.figure(figsize=(16,10))
plt.ylabel('', fontsize=12,fontname="Arial")
plt.xticks(rotation=45, fontname="Arial", fontsize=12)
plt.yticks(fontname="Arial", fontsize=12)
# sns.swarmplot(x="", y=y_ax, data=df_gr1_melt, color="blue", size=5, zorder=0)
# ax=sns.boxplot(x="", y=y_ax, data=df_gr1_melt, palette=["gray", "whitesmoke"], showfliers = False, width=[0.3], boxprops=dict(edgecolor='black',linewidth=3))
ax = sns.boxplot(x="", y="marker_signal", hue="type", data=element_pmt_melt_all, showfliers = False, palette='Pastel1')
plt.setp(ax.lines, color="grey", linewidth=0.5)
plt.setp(ax.spines.values(), color="black", linewidth=0.5)
ax = sns.swarmplot(x="", y="marker_signal", hue="type", data=element_pmt_melt_all, color='grey', size=1.5, dodge=True)

plt.savefig(r'C:\Users\libin\UCSF\MMR\enrichment/promoter_enrichment_0917.pdf', transparent=True) 
