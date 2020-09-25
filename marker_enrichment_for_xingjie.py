# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 06:42:07 2020

@author: bingkun
"""

import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib as mpl
mpl.rcParams['figure.dpi']= 300

marker_signal_all = pd.read_csv(r'C:\Users\libin\UCSF\MMR\enrichment\dataframe_melt_concat.csv', sep="\t",
                                names=["marker", "signal", "type"])

marker_list = ['ATAC','ChIP_H3K27ac','Cut_Tag_H3K27ac','CTCF','H3K4me3']
for m in marker_list:      
        dataframe_melt_concat_subset = marker_signal_all[marker_signal_all["marker"] == m]
        ax = plt.figure(figsize=(16,10))
        plt.ylabel('', fontsize=12,fontname="Arial")
        plt.xticks(rotation=45, fontname="Arial", fontsize=12)
        plt.yticks(fontname="Arial", fontsize=12)
        ax = sns.boxplot(x="marker", y="signal", hue="type", data=dataframe_melt_concat_subset, showfliers = False, palette='Pastel1')
        plt.setp(ax.lines, color="grey", linewidth=0.5)
        plt.setp(ax.spines.values(), color="black", linewidth=0.5)
        ax = sns.swarmplot(x="marker", y="signal", hue="type", data=dataframe_melt_concat_subset, color='grey', size=1.5, dodge=True)
        plt.savefig('element_enrichment_{}.pdf'.format(m), transparent=True)