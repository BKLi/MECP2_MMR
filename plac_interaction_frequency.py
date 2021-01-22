# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 22:18:57 2020

@author: bingkun
"""

import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib as mpl
mpl.rcParams['figure.dpi']= 300

#%%
element_ann = pd.read_csv(r'C:\Users\libin\UCSF\MMR\manual_annotation_complete_update_refseq_0801.csv')
element_ann = element_ann[['chr','Start','End','annotation']]
element_ann_pmt = element_ann[element_ann["annotation"] == 'promoter']

element_strong = pd.read_csv(r'C:\Users\libin\UCSF\MMR\enrichment\element_strong', sep="\t", names=["chr", "Start", "End"])
element_weak = pd.read_csv(r'C:\Users\libin\UCSF\MMR\enrichment\element_weak', sep="\t", names=["chr", "Start", "End"])
element_negative = pd.read_csv(r'C:\Users\libin\UCSF\MMR\enrichment\element_negative', sep="\t", names=["chr", "Start", "End"])

element_pos = pd.concat([element_strong, element_weak])

element_pos_pmt = pd.merge(element_ann_pmt, element_pos, on=['chr','Start','End'], how="inner")
element_pos_pmt["ID"] = ["element_pos_pmt_{}".format(i) for i in range(element_pos_pmt.shape[0])]
element_pos_pmt.to_csv(r'C:\Users\libin\UCSF\MMR\plac\element_pos_pmt',sep="\t", index=False, header=False)
element_neg_pmt = pd.merge(element_ann_pmt, element_negative, on=['chr','Start','End'], how="inner")
element_neg_pmt["ID"] = ["element_neg_pmt_{}".format(i) for i in range(element_neg_pmt.shape[0])]
element_neg_pmt.to_csv(r'C:\Users\libin\UCSF\MMR\plac\element_neg_pmt',sep="\t", index=False, header=False)

#%%
raw_full_interactions = pd.read_csv(r'C:\Users\libin\UCSF\MMR\plac\combined.5k.2.peaks', delim_whitespace=True)

full_anchor_anchor = raw_full_interactions[(raw_full_interactions["X1D_peak_bin1"] == 1) & (raw_full_interactions["X1D_peak_bin2"] == 1)]
# anchor list part one
full_anchor_anchor_anchor1 = full_anchor_anchor[["chr", "bin1_mid", "X1D_peak_bin1"]]
full_anchor_anchor_anchor2 = full_anchor_anchor[["chr", "bin2_mid", "X1D_peak_bin2"]]

full_bin1_anchor = raw_full_interactions[(raw_full_interactions["X1D_peak_bin1"] == 1) & (raw_full_interactions["X1D_peak_bin2"] == 0)]
# anchor list part two
full_bin1_anchor_anchor = full_bin1_anchor[["chr", "bin1_mid", "X1D_peak_bin1"]]

full_bin2_anchor = raw_full_interactions[(raw_full_interactions["X1D_peak_bin1"] == 0) & (raw_full_interactions["X1D_peak_bin2"] == 1)]
# anchor list part three
full_bin2_anchor_anchor = full_bin2_anchor[["chr", "bin2_mid", "X1D_peak_bin2"]]

# combine & dedup
anchor_list = pd.DataFrame(np.concatenate([full_anchor_anchor_anchor1, full_anchor_anchor_anchor2, full_bin1_anchor_anchor, full_bin2_anchor_anchor]), columns=full_anchor_anchor_anchor1.columns)
anchor_list = anchor_list.drop_duplicates()

anchor_list["start"] = anchor_list["bin1_mid"].astype('int64')
#anchor_list["start"] = anchor_list["start"].astype('int64')
anchor_list["end"] = anchor_list["bin1_mid"] + 5000
anchor_list["end"] = anchor_list["end"].astype('int64')
anchor_list = anchor_list[["chr", "start", "end"]]
anchor_list["anchor_ID"] = ["anchor_{}".format(i) for i in range(anchor_list.shape[0])]
anchor_list.to_csv(r'C:\Users\libin\UCSF\MMR\plac\anchor_list',sep="\t", index=False, header=False)

print ('number of anchors in total: {}'.format(anchor_list.shape[0]))

#%%
####### read in & process significant MAPS interaction sets

inter_sig = pd.read_csv(r'C:\Users\libin\UCSF\MMR\plac\combined.5k.2.peaks.bedpe', delim_whitespace=True)

inter_sig["inter_ID"] = ["inter_{}".format(i) for i in range(inter_sig.shape[0])]
print ("number of signif interactions: {}".format(inter_sig.shape[0]))

inter_sig_left_anchor = pd.merge(inter_sig, anchor_list, left_on=['chr1', 'start1', 'end1'], right_on = ['chr', 'start', 'end'], how="inner")
left_anchors_list = inter_sig_left_anchor[["inter_ID", 'anchor_ID']]

inter_sig_right_anchor = pd.merge(inter_sig, anchor_list, left_on=['chr2', 'start2', 'end2'], right_on = ['chr', 'start', 'end'], how="inner")
right_anchors_list = inter_sig_right_anchor[["inter_ID", 'anchor_ID']]

anchor_anchor_list = pd.merge(left_anchors_list, right_anchors_list, on="inter_ID", how="inner")
AND_interactions = inter_sig[inter_sig["inter_ID"].isin(anchor_anchor_list["inter_ID"])]
XOR_interactions = inter_sig[~inter_sig["inter_ID"].isin(anchor_anchor_list["inter_ID"])]

#%%
#reformat XOR
inter_sig_left_anchor_XOR = pd.merge(XOR_interactions, anchor_list, left_on=['chr1', 'start1', 'end1'], right_on = ['chr', 'start', 'end'], how="inner")
print ("number of signif interactions with bin1 as anchor bin: {}".format(inter_sig_left_anchor_XOR.shape[0]))

inter_sig_right_anchor_XOR = pd.merge(XOR_interactions, anchor_list, left_on=['chr2', 'start2', 'end2'], right_on = ['chr', 'start', 'end'], how="inner")
print ("number of signif interactions with bin2 as anchor bin: {}".format(inter_sig_right_anchor_XOR.shape[0]))
# rename so bin1 is always anchor bin
inter_sig_right_anchor_XOR = inter_sig_right_anchor_XOR.rename(columns={"chr1":"chr2", "start1":"start2", "end1":"end2", "chr2":"chr1", "start2":"start1", "end2":"end1"})
inter_sig_right_anchor_XOR = inter_sig_right_anchor_XOR[['chr1','start1','end1','chr2','start2','end2','count','expected','fdr','ClusterLabel','ClusterSize',
                                                         'ClusterType','ClusterNegLog10P','ClusterSummit','inter_ID','chr','start','end','anchor_ID']]

XOR_interactions_reformed = pd.concat([inter_sig_left_anchor_XOR, inter_sig_right_anchor_XOR])
XOR_interactions_reformed = XOR_interactions_reformed[['chr1','start1','end1','chr2','start2','end2','count','expected','fdr','ClusterLabel','ClusterSize',
                                                         'ClusterType','ClusterNegLog10P','ClusterSummit','inter_ID','anchor_ID']]

#%%reform AND
# collapse AND into 2, treat as 2 set of XOR interactions
inter_sig_left_anchor_AND = pd.merge(AND_interactions, anchor_list, left_on=['chr1', 'start1', 'end1'], right_on = ['chr', 'start', 'end'], how="inner")
inter_sig_right_anchor_AND = pd.merge(AND_interactions, anchor_list, left_on=['chr2', 'start2', 'end2'], right_on = ['chr', 'start', 'end'], how="inner")
inter_sig_right_anchor_AND = inter_sig_right_anchor_AND.rename(columns={"chr1":"chr2", "start1":"start2", "end1":"end2", "chr2":"chr1", "start2":"start1", "end2":"end1"})
AND_interactions_reformed = pd.concat([inter_sig_left_anchor_AND, inter_sig_right_anchor_AND])
AND_interactions_reformed = AND_interactions_reformed[['chr1','start1','end1','chr2','start2','end2','count','expected','fdr','ClusterLabel','ClusterSize',
                                                         'ClusterType','ClusterNegLog10P','ClusterSummit','inter_ID','anchor_ID']]

#%%
# anchors intersecting with control promoters or enhancer-promoters
anchor_element_neg_pmt = pd.read_csv(r'C:\Users\libin\UCSF\MMR\plac\anchor_element_neg_pmt',sep="\t", names=["chr_anc", "start_anc", "end_anc", "anchor_ID", "chr_element", "start_element", "end_elment", "type", "element_ID"])
anchor_element_pos_pmt = pd.read_csv(r'C:\Users\libin\UCSF\MMR\plac\anchor_element_pos_pmt',sep="\t", names=["chr_anc", "start_anc", "end_anc", "anchor_ID", "chr_element", "start_element", "end_elment", "type", "element_ID"])

#%%
XOR_interactions_neg_pmt = pd.merge(XOR_interactions_reformed, anchor_element_neg_pmt, on=["anchor_ID"], how="inner")
XOR_interactions_neg_pmt = XOR_interactions_neg_pmt[["inter_ID", "element_ID"]]
XOR_interactions_neg_pmt_frq = XOR_interactions_neg_pmt.groupby(["element_ID"]).size().to_frame(name = 'freq_neg').reset_index()
#AND_interactions_neg_pmt = pd.merge(AND_interactions_reformed, anchor_element_neg_pmt, on=["anchor_ID"], how="inner")

XOR_interactions_pos_pmt = pd.merge(XOR_interactions_reformed, anchor_element_pos_pmt, on=["anchor_ID"], how="inner")
XOR_interactions_pos_pmt = XOR_interactions_pos_pmt[["inter_ID", "element_ID"]]
XOR_interactions_pos_pmt_frq = XOR_interactions_pos_pmt.groupby(["element_ID"]).size().to_frame(name = 'freq_pos').reset_index()
#AND_interactions_pos_pmt = pd.merge(AND_interactions_reformed, anchor_element_pos_pmt, on=["anchor_ID"], how="inner")

#%% plot
XOR_interactions_neg_pmt_frq = XOR_interactions_neg_pmt_frq[["freq_neg"]].reset_index()["freq_neg"]
XOR_interactions_pos_pmt_frq = XOR_interactions_pos_pmt_frq[["freq_pos"]].reset_index()["freq_pos"]

freq_for_plot = pd.DataFrame()
freq_for_plot["negative"] = XOR_interactions_neg_pmt_frq
freq_for_plot["positive"] = XOR_interactions_pos_pmt_frq
freq_for_plot_melt = pd.melt(freq_for_plot, value_name="contact_frequency", var_name="type")

#%%
x = plt.figure(figsize=(10,7))
plt.ylabel('', fontsize=12,fontname="Arial")
plt.xticks(rotation=45, fontname="Arial", fontsize=12)
plt.yticks(fontname="Arial", fontsize=12)
ax = sns.swarmplot(x="type", y="contact_frequency", data=freq_for_plot_melt, color='grey', size=4, dodge=True)
#ax = sns.boxplot(x="type", y="contact_frequency", data=freq_for_plot_melt, showfliers = False, palette='Pastel1')
#plt.setp(ax.lines, color="grey", linewidth=0.5)
#plt.setp(ax.spines.values(), color="black", linewidth=0.5)
plt.savefig(r'C:\Users\libin\UCSF\MMR\plac\plac_contact_frequency_1010.pdf', transparent=True) 
