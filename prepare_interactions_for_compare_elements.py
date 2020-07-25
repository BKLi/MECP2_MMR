# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 10:48:28 2019

@author: libin
"""

import pandas as pd
import numpy as np
import seaborn as sns

WTC_interactons = pd.read_csv(r'C:\Users\libin\UCSF\MECP2\compare_element\1204_plac\WTC11_0hr_Interaction_nofilter_wtype.bedpe', delim_whitespace=True)

bin1_anchor = WTC_interactons[(WTC_interactons["X1D_peak_bin1"] == 1) & (WTC_interactons["X1D_peak_bin2"] == 0)]
bin1_ID = ["bin1_ID_{}".format(i) for i in range(bin1_anchor.shape[0])]
bin1_anchor['ID'] = bin1_ID
# extract part 1 of anchor bins
bin1_anchor_bins = bin1_anchor[["chr1", "start1", "end1", "ID"]]

bin2_anchor = WTC_interactons[(WTC_interactons["X1D_peak_bin1"] == 0) & (WTC_interactons["X1D_peak_bin2"] == 1)]
bin2_ID = ["bin2_ID_{}".format(i) for i in range(bin2_anchor.shape[0])]
bin2_anchor['ID'] = bin2_ID
## rename so bin1 is always anchor bin
# bin2_anchor_switch = bin2_anchor.rename(columns={"chr1":"chr2", "start1":"start2", "end1":"end2", "chr2":"chr1", "start2":"start1", "end2":"end1"})
bin2_anchor_switch = bin2_anchor.rename(columns={"chr1":"tmp1", "start1":"tmp2", "end1":"tmp3", "chr2":"tmp4", "start2":"tmp5", "end2":"tmp6"})
bin2_anchor_switch = bin2_anchor_switch.rename(columns={"tmp1":"chr2", "tmp2":"start2", "tmp3":"end2", "tmp4":"chr1", "tmp5":"start1", "tmp6":"end1"})
bin2_anchor_switch = bin2_anchor_switch \
[['chr1',
 'start1',
 'end1',
 'chr2',
 'start2',
 'end2', 
 'count',
 'expected',
 'fdr',
 'X1D_peak_bin1',
 'X1D_peak_bin2',
 'ID']]
# extract part 2 of anchor bins
bin2_anchor_bins = bin2_anchor_switch[["chr1", "start1", "end1", "ID"]]

# merge bin1_anchor interaction and bin2-anchor interaction
interactions_all = pd.concat([bin1_anchor, bin2_anchor_switch]) 
# check if there are duplicates
interactions_all_dedup = interactions_all.drop_duplicates(subset=['chr1','start1','end1','chr2','start2','end2'])


all_anchor_bins = pd.concat([bin1_anchor_bins, bin2_anchor_bins])
all_anchor_bins_dedup = all_anchor_bins.drop_duplicates(subset=["chr1", "start1", "end1"], keep="first")
all_anchor_bins_dedup.to_csv(r'C:\Users\libin\UCSF\MECP2\compare_element\1204_plac\WTC11_0hr_Interaction_nofilter_anchors.bed',
                             sep="\t", index=False, header=False)

# ----------------------------------------------------------------------------------------
gene_list = pd.read_csv(r'C:\Users\libin\Desktop\gene_list.bed', delim_whitespace=True, names=["chrom", "tss", "name"])
gene_list["start"] = gene_list["tss"]-1000
gene_list["end"] = gene_list["tss"]+1000
gene_list = gene_list[['chrom', 'start', 'end', 'name']]
gene_list.to_csv(r'C:\Users\libin\UCSF\MECP2\compare_element\1204_plac\gene_list.bed', sep="\t",
                 index=False, header=False)

# -----------------------------------------------------------------------------------------

peak_list = pd.read_csv(r'C:\Users\libin\UCSF\MECP2\compare_element\1204_plac\H3K4me3_peak_list.bed', 
                        delim_whitespace=True, names=["chr", "start", "end"])
peak_list.to_csv(r'C:\Users\libin\UCSF\MECP2\compare_element\1204_plac\H3K4me3_peak_list.bed', sep="\t", index=False, header=False)

# -----------------------------------------------------------------------------------------
#gene_intersect_anchors = pd.read_csv(r'C:\Users\libin\UCSF\MECP2\compare_element\1204_plac\gene_intersect_anchor.bed', sep="\t",
#                                     names=["chr1", "start1", "end1", "ID", "chr_gene", "start_gene", "end_gene", "name_gene"])

anchors_of_interest = pd.read_csv(r'C:\Users\libin\UCSF\MECP2\compare_element\1204_plac\anchors_of_interest_uniq.bed', 
                                  sep="\t", names=["chr1", "start1", "end1", "ID"])
# extract interactions whose anchor bins overlap with promoters of gene of interest
interactions_of_interest =pd.merge(interactions_all_dedup, anchors_of_interest, on=["chr1", "start1", "end1"], how="inner")
interactions_of_interest = interactions_of_interest[['chr1','start1','end1','chr2','start2','end2','count','expected','fdr']]
interactions_of_interest.to_csv(r'C:\Users\libin\UCSF\MECP2\compare_element\1204_plac\interactions_of_interest', sep="\t",
                                index=False)
interactions_of_interest_target_bins = interactions_of_interest[['chr2','start2','end2','count','expected','fdr']]
interactions_of_interest_target_bins.to_csv(r'C:\Users\libin\UCSF\MECP2\compare_element\1204_plac\interactions_of_interest_target_bins',
                                            sep="\t", index=False, header=False)




# ------------------------------- NEVERMIND --------------------------------
# check why Niko's anchor bin list doesn't match with mine
test_list = pd.read_csv(r'C:\Users\libin\UCSF\MECP2\compare_element\1204_plac\WTC11_H3K4me3_merged.q1e-4.shrt_peaks.final.bed', sep="\t",
                        names=["chr1", "start1", "end1", "col4", "col5", "col6"])

test_mine_intersect = pd.merge(test_list, all_anchor_bins_dedup, how="inner", on=["chr1", "start1", "end1"])
# LOL THE PEAKS ARE NOT BINNED IN PRIGINAL FILE - IT GOT BINNED BY RUNNING MAPS

#hg38_gtf = pd.read_csv(r'C:\Users\libin\UCSF\hg38_general\gencode.v32.annotation.gtf', delim_whitespace=True, anno)
#chr20_anchor_bins = all_anchor_bins_dedup[all_anchor_bins_dedup["chr1"] == "chr20"]











