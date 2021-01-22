# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 21:08:59 2020

@author: bingkun
"""

import pandas as pd
import numpy as np
from functools import reduce
from scipy import stats
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib as mpl
mpl.rcParams['figure.dpi']= 300

#%%
negative_homer = pd.read_csv(r'C:\Users\libin\UCSF\MMR\enrichment\homer\knownResults_negative.txt', sep="\t")
negative_homer["Motif"] = negative_homer["Motif Name"].str.extract(r'(.+?)\/')
negative_homer = negative_homer.rename(columns={'Log P-value': "negative_log_P"})
negative_homer_cut = negative_homer[['Motif', "negative_log_P"]]

strong_homer = pd.read_csv(r'C:\Users\libin\UCSF\MMR\enrichment\homer\knownResults_strong.txt', sep="\t")
strong_homer["Motif"] = strong_homer["Motif Name"].str.extract(r'(.+?)\/')
strong_homer = strong_homer.rename(columns={'Log P-value': "strong_log_P"})
strong_homer_cut = strong_homer[['Motif', "strong_log_P"]]

weak_homer = pd.read_csv(r'C:\Users\libin\UCSF\MMR\enrichment\homer\knownResults_weak.txt', sep="\t")
weak_homer["Motif"] = weak_homer["Motif Name"].str.extract(r'(.+?)\/')
weak_homer = weak_homer.rename(columns={'Log P-value': "weak_log_P"})
weak_homer_cut = weak_homer[['Motif', "weak_log_P"]]

homer_list = [negative_homer_cut, strong_homer_cut, weak_homer_cut]
homer_merged = reduce(lambda left,right: pd.merge(left,right,on=["Motif"], how='inner'), homer_list)
homer_merged["Motif"] = homer_merged["Motif"].str.upper()

#%%
expression = pd.read_csv(r'C:\Users\libin\UCSF\MECP2\RNASeq\TMM\average_expression.csv', sep=",")
expression_cut = expression[['gene_name', 'average_expression']]
expression_cut = expression_cut.drop_duplicates()
homer_merged_wExpression = pd.merge(homer_merged, expression_cut, left_on=["Motif"], right_on=["gene_name"], how="inner")
homer_merged_wExpression["diff_strong_weak"] = homer_merged_wExpression["strong_log_P"]/homer_merged_wExpression["weak_log_P"]
homer_merged_wExpression["diff_strong_neg"] = homer_merged_wExpression["strong_log_P"]/homer_merged_wExpression["negative_log_P"]
homer_merged_wExpression = homer_merged_wExpression[homer_merged_wExpression['average_expression'] != 0]

candidate_TF = homer_merged_wExpression[(homer_merged_wExpression["diff_strong_weak"] >= 2) & (homer_merged_wExpression["diff_strong_neg"] >= 2)]
candidate_TF["negative"] = "negative"
candidate_TF["strong"] = "strong"
candidate_TF["weak"] = "weak"

#%%
x = plt.figure(figsize=(10,10))
plt.yticks(rotation=90, fontname="Arial", fontsize=12)
ax = sns.heatmap(candidate_TF[['strong_log_P', 'weak_log_P', "negative_log_P"]],
                 yticklabels = candidate_TF['Motif'])
plt.savefig(r'C:\Users\libin\UCSF\MMR\enrichment\homer\homer_heatmap.pdf', transparent=True) 

#%%

fig, (ax1,ax2,ax3) = plt.subplots(ncols=3, figsize=(10,15))
plt.subplots_adjust(wspace=0)
#plt.ylabel('', fontsize=12,fontname="Arial")
#plt.xticks(fontname="Arial", fontsize=20)
#plt.yticks(fontname="Arial", fontsize=20)
#plt.tick_params(axis='both', which='major', labelsize=22)
#fig, ax = plt.subplots(figsize=(7,10))
ax1.scatter(candidate_TF["strong"], candidate_TF['Motif'], 
                     cmap="RdBu", c=candidate_TF['strong_log_P'], s=candidate_TF['average_expression']*4)
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.tick_params(labelsize=20)
               
ax2.scatter(candidate_TF['weak'], candidate_TF['Motif'], 
                     cmap="RdBu", c=candidate_TF['weak_log_P'], s=candidate_TF['average_expression']*4)
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax2.set_yticks([])
ax2.tick_params(labelsize=20)

ax3.scatter(candidate_TF['negative'], candidate_TF['Motif'], 
                     cmap="RdBu", c=candidate_TF['negative_log_P'], s=candidate_TF['average_expression']*4)
ax3.spines['right'].set_visible(False)
ax3.spines['top'].set_visible(False)
ax3.spines['left'].set_visible(False)
ax3.set_yticks([])
ax3.tick_params(labelsize=20)

plt.savefig(r'C:\Users\libin\UCSF\MMR\enrichment\homer\homer_bubble.pdf', transparent=True) 
