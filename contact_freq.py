# -*- coding: utf-8 -*-
"""
Created on Fri Oct 9 20:57:03 2020

@author: bingkun
@project: MMR
"""

import pandas as pd
import numpy as np
from scipy import stats
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib as mpl
mpl.rcParams['figure.dpi']= 300

date="1102"
#%% 
chrom_list = ["chr3", "chr20", "chrX", "chr7"]

element_weak_bin_list = []
element_strong_bin_list = []
promoter_bin_list = []

negative_freq_list = []
weak_freq_list = []
strong_freq_list = []

negative_sumCount_list = []
weak_sumCount_list = []
strong_sumCount_list = []

negative_element_promoter_dist_list = []
weak_element_promoter_dist_list = []
strong_element_promoter_dist_list = []

for chrom in chrom_list:
        #chrom="chr3"
        # group1: element in bin1, promoter in bin2
        promter_bin2 = pd.read_csv(r'C:\Users\libin\UCSF\MMR\hic\analysis\intersect\promoter_bin2_{}'.format(chrom),
                                   sep="\t", names=["inter_chr", "inter_start2", "inter_end2", "inter_ID", "count", "pmt_chr", "pmt_start", "pmt_end", "tx_id", "gene_name"])
        promter_bin2["count"] = promter_bin2["count"].round(1)
        promoter_bin_list.append(promter_bin2)
        
        element_negative_bin1 = pd.read_csv(r'C:\Users\libin\UCSF\MMR\hic\analysis\intersect\element_negative_bin1_{}'.format(chrom),
                                            sep="\t", names=["inter_chr", "inter_start1", "inter_end1", "inter_ID", "count", "element_chr", "element_start", "element_end"])
        element_negative_bin1["count"] = element_negative_bin1["count"].round(1)
        negative_element_bin1_promoter_bin2 = pd.merge(element_negative_bin1, promter_bin2, on=["inter_chr", "inter_ID", "count"], how="inner")
        # remove items with same interaction, same gene => element contacting different tx of same gene by same interactions is considered caontacting the gene only once
        negative_element_bin1_promoter_bin2 = negative_element_bin1_promoter_bin2.drop_duplicates(subset=["inter_ID", "gene_name", "element_chr", "element_start", "element_end"], keep="first")
        negative_element_bin1_promoter_bin2_inter = negative_element_bin1_promoter_bin2[["inter_chr", "inter_start1", "inter_end1", "inter_ID", "inter_start2", "inter_end2"]]
        negative_element_bin1_promoter_bin2 = negative_element_bin1_promoter_bin2[["inter_ID", "count", "gene_name", "element_chr", "element_start", "element_end"]]
        # group2: element in bin2, promoter in bin1
        promter_bin1 = pd.read_csv(r'C:\Users\libin\UCSF\MMR\hic\analysis\intersect\promoter_bin1_{}'.format(chrom),
                                   sep="\t", names=["inter_chr", "inter_start1", "inter_end1", "inter_ID", "count", "pmt_chr", "pmt_start", "pmt_end", "tx_id", "gene_name"])
        promoter_bin_list.append(promter_bin1)        
        promter_bin1["count"] = promter_bin1["count"].round(1)        
        element_negative_bin2 = pd.read_csv(r'C:\Users\libin\UCSF\MMR\hic\analysis\intersect\element_negative_bin2_{}'.format(chrom),
                                            sep="\t", names=["inter_chr", "inter_start2", "inter_end2", "inter_ID", "count", "element_chr", "element_start", "element_end"])
        element_negative_bin2["count"] = element_negative_bin2["count"].round(1)                
        negative_element_bin2_promoter_bin1 = pd.merge(element_negative_bin2, promter_bin1, on=["inter_chr", "inter_ID", "count"], how="inner")
        negative_element_bin2_promoter_bin1 = negative_element_bin2_promoter_bin1.drop_duplicates(subset=["inter_ID", "gene_name", "element_chr", "element_start", "element_end"], keep="first")
        negative_element_bin2_promoter_bin1_inter = negative_element_bin2_promoter_bin1[["inter_chr", "inter_start1", "inter_end1", "inter_ID", "inter_start2", "inter_end2"]]
        negative_element_bin2_promoter_bin1 = negative_element_bin2_promoter_bin1[["inter_ID", "count", "gene_name", "element_chr", "element_start", "element_end"]]
        
        #negative_element_bin1_promoter_bin2 = negative_element_bin1_promoter_bin2.dropna()
        #negative_element_bin2_promoter_bin1 = negative_element_bin2_promoter_bin1.dropna()
        
        negative_element_promoter_concat = pd.concat([negative_element_bin1_promoter_bin2, negative_element_bin2_promoter_bin1])
        negative_element_promoter_concat = negative_element_promoter_concat.drop_duplicates()
        print(negative_element_promoter_concat.shape[0])
        
        negative_element_promoter_inter_concat = pd.concat([negative_element_bin1_promoter_bin2_inter,negative_element_bin2_promoter_bin1_inter])
        negative_element_promoter_inter_concat_dist = pd.merge(negative_element_promoter_inter_concat, negative_element_promoter_concat, on=["inter_ID"], how="inner")
        negative_element_promoter_inter_concat_dist = negative_element_promoter_inter_concat_dist.drop_duplicates()
        print(negative_element_promoter_inter_concat_dist.shape[0])
        
        negative_element_promoter_dist_list.append(negative_element_promoter_inter_concat_dist)
        # count occurance of each element : https://stackoverflow.com/questions/38933071/group-by-two-columns-and-count-the-occurrences-of-each-combination-in-pandas
        #contact_freq_count_neg = negative_element_promoter_concat.groupby(["element_chr", "element_start", "element_end"]).size().to_frame(name = 'freq_neg').reset_index()
        #contact_freq_count_neg["element_ID"] = ["{}_{}".format(chrom, i) for i in range(contact_freq_count_neg.shape[0])]
        #negative_freq_list.append(contact_freq_count_neg)
        
        # take sum of count for each element -- **interactivity**
        negative_element_promoter_concat_sumCount = negative_element_promoter_concat.groupby(["element_chr", "element_start", "element_end"]).agg(sum).reset_index()
        negative_element_promoter_concat_sumCount["element_ID"] = ["{}_{}".format(chrom, i) for i in range(negative_element_promoter_concat_sumCount.shape[0])]
        negative_sumCount_list.append(negative_element_promoter_concat_sumCount)
        
        
        element_weak_bin1 = pd.read_csv(r'C:\Users\libin\UCSF\MMR\hic\analysis\intersect\element_weak_bin1_{}'.format(chrom),
                                            sep="\t", names=["inter_chr", "inter_start1", "inter_end1", "inter_ID", "count", "element_chr", "element_start", "element_end"])
        element_weak_bin1["count"] = element_weak_bin1["count"].round(1)                
        weak_element_bin1_promoter_bin2 = pd.merge(element_weak_bin1, promter_bin2, on=["inter_chr", "inter_ID", "count"], how="inner")
        # remove items with same interaction, same gene => element contacting different tx of same gene by same interactions is considered caontacting the gene only once
        weak_element_bin1_promoter_bin2 = weak_element_bin1_promoter_bin2.drop_duplicates(subset=["inter_ID", "gene_name", "element_chr", "element_start", "element_end"], keep="first")
        weak_element_bin1_promoter_bin2_inter = weak_element_bin1_promoter_bin2[["inter_chr", "inter_start1", "inter_end1", "inter_ID", "inter_start2", "inter_end2"]]
        weak_element_bin1_promoter_bin2 = weak_element_bin1_promoter_bin2[["inter_ID", "count", "gene_name", "element_chr", "element_start", "element_end"]]
        # group2: element in bin2, promoter in bin1
        element_weak_bin2 = pd.read_csv(r'C:\Users\libin\UCSF\MMR\hic\analysis\intersect\element_weak_bin2_{}'.format(chrom),
                                            sep="\t", names=["inter_chr", "inter_start2", "inter_end2", "inter_ID", "count", "element_chr", "element_start", "element_end"])
        element_weak_bin2["count"] = element_weak_bin2["count"].round(1)                
       
        weak_element_bin2_promoter_bin1 = pd.merge(element_weak_bin2, promter_bin1, on=["inter_chr", "inter_ID", "count"], how="inner")
        weak_element_bin2_promoter_bin1 = weak_element_bin2_promoter_bin1.drop_duplicates(subset=["inter_ID", "gene_name", "element_chr", "element_start", "element_end"], keep="first")
        weak_element_bin2_promoter_bin1_inter = weak_element_bin2_promoter_bin1[["inter_chr", "inter_start1", "inter_end1", "inter_ID", "inter_start2", "inter_end2"]]
        weak_element_bin2_promoter_bin1 = weak_element_bin2_promoter_bin1[["inter_ID", "count", "gene_name", "element_chr", "element_start", "element_end"]]
        
        #weak_element_bin1_promoter_bin2 = weak_element_bin1_promoter_bin2.dropna()
        #weak_element_bin2_promoter_bin1 = weak_element_bin2_promoter_bin1.dropna()
        
        weak_element_promoter_concat = pd.concat([weak_element_bin1_promoter_bin2, weak_element_bin2_promoter_bin1])
        weak_element_promoter_concat = weak_element_promoter_concat.drop_duplicates()
        
        weak_element_promoter_inter_concat = pd.concat([weak_element_bin1_promoter_bin2_inter,weak_element_bin2_promoter_bin1_inter])
        weak_element_promoter_inter_concat_dist = pd.merge(weak_element_promoter_inter_concat, weak_element_promoter_concat, on=["inter_ID"], how="inner")        
        weak_element_promoter_inter_concat_dist = weak_element_promoter_inter_concat_dist.drop_duplicates()
        weak_element_promoter_dist_list.append(weak_element_promoter_inter_concat_dist)
        
        element_weak_bin_list.append(element_weak_bin1)
        element_weak_bin_list.append(element_weak_bin2)
        # count occurance of each element : https://stackoverflow.com/questions/38933071/group-by-two-columns-and-count-the-occurrences-of-each-combination-in-pandas
        #contact_freq_count_weak = weak_element_promoter_concat.groupby(["element_chr", "element_start", "element_end"]).size().to_frame(name = 'freq_weak').reset_index()
        #contact_freq_count_weak["element_ID"] = ["{}_{}".format(chrom, i) for i in range(contact_freq_count_weak.shape[0])]
        #weak_freq_list.append(contact_freq_count_weak)       
        
        weak_element_promoter_concat_sumCount = weak_element_promoter_concat.groupby(["element_chr", "element_start", "element_end"]).agg(sum).reset_index()
        weak_element_promoter_concat_sumCount["element_ID"] = ["{}_{}".format(chrom, i) for i in range(weak_element_promoter_concat_sumCount.shape[0])]
        weak_sumCount_list.append(weak_element_promoter_concat_sumCount)
        
# --------------------------------------------------------------------------------------------------------        
        element_strong_bin1 = pd.read_csv(r'C:\Users\libin\UCSF\MMR\hic\analysis\intersect\element_strong_bin1_{}'.format(chrom),
                                            sep="\t", names=["inter_chr", "inter_start1", "inter_end1", "inter_ID", "count", "element_chr", "element_start", "element_end"])
        element_strong_bin1["count"] = element_strong_bin1["count"].round(1)                

        strong_element_bin1_promoter_bin2 = pd.merge(element_strong_bin1, promter_bin2, on=["inter_chr", "inter_ID", "count"], how="inner")
        # remove items with same interaction, same gene => element contacting different tx of same gene by same interactions is considered caontacting the gene only once
        strong_element_bin1_promoter_bin2 = strong_element_bin1_promoter_bin2.drop_duplicates(subset=["inter_ID", "gene_name", "element_chr", "element_start", "element_end"], keep="first")
        strong_element_bin1_promoter_bin2_inter = strong_element_bin1_promoter_bin2[["inter_chr", "inter_start1", "inter_end1", "inter_ID", "inter_start2", "inter_end2"]]
        strong_element_bin1_promoter_bin2 = strong_element_bin1_promoter_bin2[["inter_ID", "count", "gene_name", "element_chr", "element_start", "element_end"]]
        # group2: element in bin2, promoter in bin1
        element_strong_bin2 = pd.read_csv(r'C:\Users\libin\UCSF\MMR\hic\analysis\intersect\element_strong_bin2_{}'.format(chrom),
                                            sep="\t", names=["inter_chr", "inter_start2", "inter_end2", "inter_ID", "count", "element_chr", "element_start", "element_end"])
        element_strong_bin2["count"] = element_strong_bin2["count"].round(1)                
        strong_element_bin2_promoter_bin1 = pd.merge(element_strong_bin2, promter_bin1, on=["inter_chr", "inter_ID", "count"], how="inner")
        strong_element_bin2_promoter_bin1 = strong_element_bin2_promoter_bin1.drop_duplicates(subset=["inter_ID", "gene_name", "element_chr", "element_start", "element_end"], keep="first")
        strong_element_bin2_promoter_bin1_inter = strong_element_bin2_promoter_bin1[["inter_chr", "inter_start1", "inter_end1", "inter_ID", "inter_start2", "inter_end2"]]
        strong_element_bin2_promoter_bin1 = strong_element_bin2_promoter_bin1[["inter_ID", "count", "gene_name", "element_chr", "element_start", "element_end"]]
        
        #strong_element_bin1_promoter_bin2 = strong_element_bin1_promoter_bin2.dropna()
        #strong_element_bin2_promoter_bin1 = strong_element_bin2_promoter_bin1.dropna()
        
        strong_element_promoter_concat = pd.concat([strong_element_bin1_promoter_bin2, strong_element_bin2_promoter_bin1])
        strong_element_promoter_concat = strong_element_promoter_concat.drop_duplicates()
        
        strong_element_promoter_inter_concat = pd.concat([strong_element_bin1_promoter_bin2_inter,strong_element_bin2_promoter_bin1_inter])
        strong_element_promoter_inter_concat_dist = pd.merge(strong_element_promoter_inter_concat, strong_element_promoter_concat, on=["inter_ID"], how="inner")        
        strong_element_promoter_inter_concat_dist = strong_element_promoter_inter_concat_dist.drop_duplicates()
        strong_element_promoter_dist_list.append(strong_element_promoter_inter_concat_dist)
        
        # count occurance of each element : https://stackoverflow.com/questions/38933071/group-by-two-columns-and-count-the-occurrences-of-each-combination-in-pandas
        #contact_freq_count_strong = strong_element_promoter_concat.groupby(["element_chr", "element_start", "element_end"]).size().to_frame(name = 'freq_strong').reset_index()
        #contact_freq_count_strong["element_ID"] = ["{}_{}".format(chrom, i) for i in range(contact_freq_count_strong.shape[0])]
        #strong_freq_list.append(contact_freq_count_strong)
        element_strong_bin_list.append(element_strong_bin1)
        element_strong_bin_list.append(element_strong_bin2)
        
        strong_element_promoter_concat_sumCount = strong_element_promoter_concat.groupby(["element_chr", "element_start", "element_end"]).agg(sum).reset_index()
        strong_element_promoter_concat_sumCount["element_ID"] = ["{}_{}".format(chrom, i) for i in range(strong_element_promoter_concat_sumCount.shape[0])]
        strong_sumCount_list.append(strong_element_promoter_concat_sumCount)
 
# -----------------------------------------------------------------------------
promoter_bin = pd.concat(promoter_bin_list)

promoter_bin_grp1 = promoter_bin.dropna(subset=["inter_start1"])
promoter_bin_grp1 = promoter_bin_grp1.drop_duplicates(subset=["pmt_chr", "pmt_start", "pmt_end", "inter_start1", "inter_end1"])
promoter_bin_grp1 = promoter_bin_grp1.rename(columns= {"inter_start1": "inter_start", "inter_end1": "inter_end"})
promoter_bin_grp1 = promoter_bin_grp1[["pmt_chr", "pmt_start", "pmt_end", "inter_start", "inter_end"]]

promoter_bin_grp2 = promoter_bin.dropna(subset=["inter_start2"])
promoter_bin_grp2 = promoter_bin_grp2.drop_duplicates(subset=["pmt_chr", "pmt_start", "pmt_end", "inter_start2", "inter_end2"])
promoter_bin_grp2 = promoter_bin_grp2.rename(columns= {"inter_start2": "inter_start", "inter_end2": "inter_end"})
promoter_bin_grp2 = promoter_bin_grp2[["pmt_chr", "pmt_start", "pmt_end", "inter_start", "inter_end"]]

promoter_bin_update = pd.concat([promoter_bin_grp1, promoter_bin_grp2])
promoter_bin_update = promoter_bin_update.drop_duplicates()

promoter_bin_update["inter_start"] = promoter_bin_update["inter_start"].apply(int)
promoter_bin_update["inter_end"] = promoter_bin_update["inter_end"].apply(int)

# -----------------------------------------------------------------------------        

element_weak_bin = pd.concat(element_weak_bin_list) 
element_weak_bin_grp1 = element_weak_bin.dropna(subset=["inter_start1"])
element_weak_bin_grp1 = element_weak_bin_grp1.drop_duplicates(subset=["element_chr", "element_start", "element_end", "inter_start1", "inter_end1"])
element_weak_bin_grp1 = element_weak_bin_grp1.rename(columns= {"inter_start1": "inter_start", "inter_end1": "inter_end"})
element_weak_bin_grp1 = element_weak_bin_grp1[["element_chr", "element_start", "element_end", "inter_start", "inter_end"]]

element_weak_bin_grp2 = element_weak_bin.dropna(subset=["inter_start2"])
element_weak_bin_grp2 = element_weak_bin_grp2.drop_duplicates(subset=["element_chr", "element_start", "element_end", "inter_start1", "inter_end1"])
element_weak_bin_grp2 = element_weak_bin_grp2.rename(columns= {"inter_start2": "inter_start", "inter_end2": "inter_end"})
element_weak_bin_grp2 = element_weak_bin_grp2[["element_chr", "element_start", "element_end", "inter_start", "inter_end"]]

element_weak_bin_update = pd.concat([element_weak_bin_grp1, element_weak_bin_grp2])
element_weak_bin_update = element_weak_bin_update.drop_duplicates()
element_weak_bin_update["inter_start"] = element_weak_bin_update["inter_start"].apply(int)
element_weak_bin_update["inter_end"] = element_weak_bin_update["inter_end"].apply(int)

weak_concat = pd.concat(weak_sumCount_list)
element_weak = pd.read_csv(r'C:\Users\libin\UCSF\MMR\enrichment\element_weak', sep="\t", names=["chr", "start", "end"])
element_weak_missing = pd.merge(weak_concat, element_weak, left_on = ['element_chr', 'element_start', 'element_end'], right_on=["chr", "start", "end"], indicator=True, how="outer")
element_weak_missing = element_weak_missing[element_weak_missing["chr"] != "chr2"] 
element_weak_missing = element_weak_missing[element_weak_missing["_merge"] == "right_only"]

element_weak_missing_cut = element_weak_missing[["chr", "start", "end"]]
element_weak_missing_bin = pd.merge(element_weak_bin_update, element_weak_missing_cut, left_on = ['element_chr', 'element_start', 'element_end'], right_on=["chr", "start", "end"], how="inner")

element_weak_missing_bin_promoter_bin = pd.merge(element_weak_missing_bin, promoter_bin_update, left_on=["element_chr"], right_on=["pmt_chr"], how="inner")
element_weak_missing_bin_promoter_bin["dist"] = abs(element_weak_missing_bin_promoter_bin["inter_start_y"] - element_weak_missing_bin_promoter_bin["inter_start_x"])
element_weak_missing_bin_promoter_bin = element_weak_missing_bin_promoter_bin[['element_chr', 'element_start', 'element_end', 'dist']]
element_weak_missing_bin_promoter_bin = element_weak_missing_bin_promoter_bin.groupby(['element_chr', 'element_start', 'element_end']).agg(np.mean).reset_index()
element_weak_missing_bin_promoter_bin["type"] = "weak"
element_weak_missing_bin_promoter_bin["count"] = 0

#negative_count_concat = pd.concat(negative_sumCount_list)[["count"]].reset_index()["count"]
weak_count_concat = pd.concat(weak_sumCount_list)[["count"]].reset_index()["count"]
weak_count_concat = weak_count_concat.append(pd.Series([0]*element_weak_missing.shape[0]))
weak_count_concat = weak_count_concat.reset_index()[0]


#------------------------------------------------------------------------------


element_strong_bin = pd.concat(element_strong_bin_list) 
element_strong_bin_grp1 = element_strong_bin.dropna(subset=["inter_start1"])
element_strong_bin_grp1 = element_strong_bin_grp1.drop_duplicates(subset=["element_chr", "element_start", "element_end", "inter_start1", "inter_end1"])
element_strong_bin_grp1 = element_strong_bin_grp1.rename(columns= {"inter_start1": "inter_start", "inter_end1": "inter_end"})
element_strong_bin_grp1 = element_strong_bin_grp1[["element_chr", "element_start", "element_end", "inter_start", "inter_end"]]

element_strong_bin_grp2 = element_strong_bin.dropna(subset=["inter_start2"])
element_strong_bin_grp2 = element_strong_bin_grp2.drop_duplicates(subset=["element_chr", "element_start", "element_end", "inter_start1", "inter_end1"])
element_strong_bin_grp2 = element_strong_bin_grp2.rename(columns= {"inter_start2": "inter_start", "inter_end2": "inter_end"})
element_strong_bin_grp2 = element_strong_bin_grp2[["element_chr", "element_start", "element_end", "inter_start", "inter_end"]]

element_strong_bin_update = pd.concat([element_strong_bin_grp1, element_strong_bin_grp2])
element_strong_bin_update = element_strong_bin_update.drop_duplicates()
element_strong_bin_update["inter_start"] = element_strong_bin_update["inter_start"].apply(int)
element_strong_bin_update["inter_end"] = element_strong_bin_update["inter_end"].apply(int)

strong_concat = pd.concat(strong_sumCount_list)
element_strong = pd.read_csv(r'C:\Users\libin\UCSF\MMR\enrichment\element_strong', sep="\t", names=["chr", "start", "end"])
element_strong_missing = pd.merge(strong_concat, element_strong, left_on = ['element_chr', 'element_start', 'element_end'], right_on=["chr", "start", "end"], indicator=True, how="outer")
element_strong_missing = element_strong_missing[element_strong_missing["chr"] != "chr2"] 
element_strong_missing = element_strong_missing[element_strong_missing["_merge"] == "right_only"]

element_strong_missing_cut = element_strong_missing[["chr", "start", "end"]]
element_strong_missing_bin = pd.merge(element_strong_bin_update, element_strong_missing_cut, left_on = ['element_chr', 'element_start', 'element_end'], right_on=["chr", "start", "end"], how="inner")

element_strong_missing_bin_promoter_bin = pd.merge(element_strong_missing_bin, promoter_bin_update, left_on=["element_chr"], right_on=["pmt_chr"], how="inner")
element_strong_missing_bin_promoter_bin["dist"] = abs(element_strong_missing_bin_promoter_bin["inter_start_y"] - element_strong_missing_bin_promoter_bin["inter_start_x"])
element_strong_missing_bin_promoter_bin = element_strong_missing_bin_promoter_bin[['element_chr', 'element_start', 'element_end', 'dist']]
element_strong_missing_bin_promoter_bin = element_strong_missing_bin_promoter_bin.groupby(['element_chr', 'element_start', 'element_end']).agg(np.mean).reset_index()
element_strong_missing_bin_promoter_bin["type"] = "strong"
element_strong_missing_bin_promoter_bin["count"] = 0
          
strong_count_concat = pd.concat(strong_sumCount_list)[["count"]].reset_index()["count"]
strong_count_concat = strong_count_concat.append(pd.Series([0]*element_strong_missing.shape[0]))
strong_count_concat = strong_count_concat.reset_index()[0]

#%%
count_for_plot = pd.DataFrame()
#count_for_plot["negative"] = negative_count_concat
count_for_plot["weak"] = weak_count_concat
count_for_plot["strong"] = strong_count_concat
count_for_plot_melt = pd.melt(count_for_plot, value_name="sum_count", var_name="type")
weak_strong_count = stats.kruskal(weak_count_concat, strong_count_concat, nan_policy="omit")[1]

#%%
x = plt.figure(figsize=(16,10))
plt.ylabel('', fontsize=12,fontname="Arial")
plt.xticks(rotation=45, fontname="Arial", fontsize=12)
plt.yticks(fontname="Arial", fontsize=12)
ax = sns.swarmplot(x="type", y="sum_count", data=count_for_plot_melt, color='grey', size=4, dodge=True)
ax = sns.boxplot(x="type", y="sum_count", data=count_for_plot_melt, showfliers = False, palette='Pastel1')
plt.setp(ax.lines, color="grey", linewidth=0.5)
plt.setp(ax.spines.values(), color="black", linewidth=0.5)
plt.savefig(r'C:\Users\libin\UCSF\MMR\hic\analysis\hic_contact_freq_{}.pdf'.format(date), transparent=True) 

#%%

element_strong_promoter_dist_all = pd.concat(strong_element_promoter_dist_list,sort=False)
element_strong_promoter_dist_all["dist"] = element_strong_promoter_dist_all["inter_start2"] - element_strong_promoter_dist_all["inter_start1"]
element_strong_promoter_dist_all["type"] = "strong"
#element_strong_promoter_dist_all_freq = element_strong_promoter_dist_all.groupby(["element_chr", "element_start", "element_end"]).size().to_frame(name = 'freq_strong').reset_index()

element_weak_promoter_dist_all = pd.concat(weak_element_promoter_dist_list,sort=False)
element_weak_promoter_dist_all["dist"] = element_weak_promoter_dist_all["inter_start2"] - element_weak_promoter_dist_all["inter_start1"]
element_weak_promoter_dist_all["type"] = "weak"

element_negative_promoter_dist_all = pd.concat(negative_element_promoter_dist_list,sort=False)
element_negative_promoter_dist_all["dist"] = element_negative_promoter_dist_all["inter_start2"] - element_negative_promoter_dist_all["inter_start1"]
element_negative_promoter_dist_all["type"] = "negative"

element_promoter_all = pd.concat([element_strong_promoter_dist_all, element_weak_promoter_dist_all], sort=False)
#element_promoter_all_frq = element_promoter_all.groupby(["element_chr", "element_start", "element_end"]).size().to_frame(name = 'freq_all').reset_index()
element_promoter_all = element_promoter_all[['element_chr', 'element_start', 'element_end', 'dist', 'type', 'count']]
#
#%%
#element_promoter_all_freq = element_promoter_all.groupby(["element_chr", "element_start", "element_end", "dist"]).size().to_frame(name = 'freq').reset_index()
#element_promoter_all_freq_annot = pd.merge(element_promoter_all_freq, element_promoter_all, on=["element_chr", "element_start", "element_end", "dist"], how="inner")
#element_promoter_all_freq_annot = element_promoter_all_freq_annot.drop_duplicates(subset=["element_chr", "element_start", "element_end", "dist"])
#element_promoter_all_freq_annot = element_promoter_all_freq_annot[["element_chr", "element_start", "element_end", "dist", "type", "freq"]]
#element_promoter_all_freq_annot_positive = element_promoter_all_freq_annot[element_promoter_all_freq_annot["type"] != "negative"]


element_promoter_all_count = element_promoter_all.groupby(["element_chr", "element_start", "element_end"]).agg({"count":"sum", "dist":np.mean}).reset_index()
element_promoter_all_count_annot = pd.merge(element_promoter_all_count, element_promoter_all, on=["element_chr", "element_start", "element_end"], how="inner")
element_promoter_all_count_annot = element_promoter_all_count_annot[["element_chr", "element_start", "element_end", "dist_x", "type", "count_x"]]
element_promoter_all_count_annot = element_promoter_all_count_annot.drop_duplicates(subset=["element_chr", "element_start", "element_end"])
#element_promoter_all_count_annot_positive = element_promoter_all_count_annot[element_promoter_all_count_annot["type"] != "negative"]
element_promoter_all_count_annot = element_promoter_all_count_annot.rename(columns={"count_x":"count", "dist_x":"dist"})
element_promoter_all_count_annot_update = pd.concat([element_promoter_all_count_annot, element_weak_missing_bin_promoter_bin, element_strong_missing_bin_promoter_bin])
element_promoter_all_count_annot_update.to_csv(r'C:\Users\libin\UCSF\MMR\hic\analysis\element_contact_freq_HiC.csv', sep=",", index=False)

element_promoter_all_count_annot_update["log_dist"] = element_promoter_all_count_annot_update['dist'].apply(np.log2)

#%%
x = plt.figure(figsize=(10,7))
plt.ylabel('', fontsize=12,fontname="Arial")
plt.xticks(rotation=45, fontname="Arial", fontsize=12)
plt.yticks(fontname="Arial", fontsize=12)
#ax = sns.boxplot(x="dist_bin", y="freq", hue="type", data=element_promoter_all_grouped_frq_combined, palette="Pastel2")
ax = sns.lmplot(x="log_dist", y="count", hue="type", data=element_promoter_all_count_annot_update)

plt.savefig(r'C:\Users\libin\UCSF\MMR\hic\analysis\hic_contact_freq_withDistance_{}.pdf'.format(date), transparent=True) 

#%%
### group by distance and count occurance in each distance group

dist_cut = 150000
dist_range = np.arange(0, element_promoter_all["dist"].max()+dist_cut, dist_cut)
dist_name = np.arange(0, len(dist_range)-1, 1)

element_promoter_all["dist_bin"] = pd.cut(element_promoter_all["dist"], bins=dist_range, labels=dist_name)


element_promoter_all_grouped_frq = []
element_promoter_all_grouped_count = []

element_promoter_all_grouped = element_promoter_all.groupby(['dist_bin'])
for name, group in element_promoter_all_grouped:
        #print(group)
        #group_freq = group.groupby(["element_chr", "element_start", "element_end"]).size().to_frame(name = 'freq').reset_index()
        #print(group_freq)
        #group_freq_annot = pd.merge(group_freq, group, on=["element_chr", "element_start", "element_end"], how="inner")
        #print(group_freq_annot)
        #element_promoter_all_grouped_frq.append(group_freq_annot)
        
        group_count = group.groupby(["element_chr", "element_start", "element_end"]).agg(sum).reset_index()
        group_count_annot = pd.merge(group_count, group, on=["element_chr", "element_start", "element_end"], how="inner")
        element_promoter_all_grouped_count.append(group_count_annot)
        

element_promoter_all_grouped_count_combined = pd.concat(element_promoter_all_grouped_count)
element_promoter_all_grouped_count_combined = element_promoter_all_grouped_count_combined[element_promoter_all_grouped_count_combined["type"] != "negative"]

element_promoter_all_grouped_count_combined_mean = element_promoter_all_grouped_count_combined.groupby(["element_chr", "element_start", "element_end", "type"]).agg({'dist': np.mean, 'count': 'sum'}).reset_index()

element_promoter_all_grouped_count_combined = element_promoter_all_grouped_count_combined.drop_duplicates(subset=["element_chr", "element_start", "element_end", "dist_bin"])
element_promoter_all_grouped_count_combined = element_promoter_all_grouped_count_combined[["element_chr", "element_start", "element_end", "count", "dist", "dist_bin", "type"]]

#%%
x = plt.figure(figsize=(16,12))
plt.ylabel('', fontsize=12,fontname="Arial")
plt.xticks(rotation=45, fontname="Arial", fontsize=12)
plt.yticks(fontname="Arial", fontsize=12)
#ax = sns.boxplot(x="dist_bin", y="freq", hue="type", data=element_promoter_all_grouped_frq_combined, palette="Pastel2")
ax = sns.scatterplot(x="count", y="dist", hue="type", data=element_promoter_all_grouped_count_combined_mean)
plt.savefig(r'C:\Users\libin\UCSF\MMR\hic\analysis\hic_count_contact_mean_distance.pdf'.format(dist_cut), transparent=True) 

#%%
x = plt.figure(figsize=(16,12))
plt.ylabel('', fontsize=12,fontname="Arial")
plt.xticks(rotation=45, fontname="Arial", fontsize=12)
plt.yticks(fontname="Arial", fontsize=12)
#ax = sns.boxplot(x="dist_bin", y="count", hue="type", data=element_promoter_all_grouped_count_combined, palette="Pastel2")
ax = sns.swarmplot(x="dist_bin", y="count", hue="type", data=element_promoter_all_grouped_count_combined, size=4, dodge=True)
plt.setp(ax.lines, color="grey", linewidth=0.5)
plt.setp(ax.spines.values(), color="black", linewidth=0.5)
plt.savefig(r'C:\Users\libin\UCSF\MMR\hic\analysis\hic_count_contact_distanceBinned_{}.pdf'.format(dist_cut), transparent=True) 

#%%
x = plt.figure(figsize=(16,12))
plt.ylabel('', fontsize=12,fontname="Arial")
plt.xticks(rotation=45, fontname="Arial", fontsize=12)
plt.yticks(fontname="Arial", fontsize=12)
ax = sns.lmplot(x="dist_bin", y="count", hue="type", data=element_promoter_all_grouped_count_combined)
plt.savefig(r'C:\Users\libin\UCSF\MMR\hic\analysis\hic_count_contact_reg_distanceBinned_{}.pdf'.format(dist_cut), transparent=True) 

#%%
x = plt.figure(figsize=(5,5))
element_promoter_all_positive = element_promoter_all[element_promoter_all["type"] != "negative"]
sns.displot(data=element_promoter_all_positive, x="dist", binwidth=10000, hue="type")
plt.savefig(r'C:\Users\libin\UCSF\MMR\hic\analysis\hic_count_distanceDistribution.pdf', transparent=True) 
