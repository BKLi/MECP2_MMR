# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 20:00:09 2020

@author: bingkun
@project : MMR

Assign genomic annotations to MMR element (ATAC-Seq peaks)

        
"""

#%%
 
#previous steps:
#       1) add ID to element
#       2)
#         if chipseeker:
#                chipseeker_MMR.R
#         if bedtools:
#                bedtools intersect -a mmr_elements_withID -b WTC11_H3K4me3_merged.q1e-4.shrt_peaks.final.bed -wa -wb > MMR_element.H4K4me3.intersect
#                bedtools intersect -a mmr_elements_withID -b hg38_exon_by_tx -wa -wb > MMR_element.exon.intersect
#                bedtools intersect -a mmr_elements_withID -b hg38_intron_by_tx -wa -wb > MMR_element.intron.intersect
#                bedtools intersect -a gencode.v32.byTX.500bp.promoter -b WTC11_H3K4me3_merged.q1e-4.shrt_peaks.final.bed -wa -wb > promoter.H3K4me3.intersect
                                  

#%%
import sys
import pandas as pd
import numpy as np
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt
from functools import reduce
import itertools
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

#%%
chipseeker = False # if annotation is chipseeker-based
bedtools = True # if annotation is bedtools-intersect based
source = 'RefSeq'
date=''

#%%
element = pd.read_csv(r'C:\Users\libin\UCSF\MMR\mmr_elements_withID', sep="\t")

hg38_geneID2Name = pd.read_csv(r'C:\Users\libin\UCSF\hg38_general\gencode.v32.geneID2geneNameType', sep="\t").dropna(axis=0)
hg38_tx2gene = pd.read_csv(r'C:\Users\libin\UCSF\hg38_general\gencode.v32.transcript2gene', sep="\t").dropna(axis=0)

#%%
if bedtools:

        promoter_intersect_H3K4me3_infile = r'C:\Users\libin\UCSF\MMR\promoter.H3K4me3.intersect'
        H3K4me3_intersect_infile = r'C:\Users\libin\UCSF\MMR\MMR_element.H4K4me3.intersect'
        intron_intersect_infile = r'C:\Users\libin\UCSF\MMR\MMR_element.intron.intersect'
        exon_intersect_infile = r'C:\Users\libin\UCSF\MMR\MMR_element.exon.intersect'
        
        # H3K4me3 linked to promoter regions they intersect with
        promoter_intersect_H3K4me3 = pd.read_csv(promoter_intersect_H3K4me3_infile, sep="\t",
                                                 names=["pmt_chr", "pmt_start", "pmt_end", "tx_id", "gene_id", "gene_name", "peak_chr", "peak_start", "peak_end", "strand", "peak_ID", "depr"])
        #de-dup : item dropped if H3K4me3 peak intersect with same gene but dif tx
        promoter_intersect_H3K4me3 = promoter_intersect_H3K4me3.drop_duplicates(subset=["gene_id", "gene_name", "peak_chr", "peak_start", "peak_end", "strand", "peak_ID", "depr"])
        # group by H3K4me3 peak (each peak may cover more than one promoter reions)
        promoter_intersect_H3K4me3 = promoter_intersect_H3K4me3[["gene_name", "peak_ID"]]
        promoter_intersect_H3K4me3 = promoter_intersect_H3K4me3. groupby(['peak_ID'], as_index=False).agg(",".join).reset_index().drop(["index"], axis=1)
        
        
        H3K4me3_intersect = pd.read_csv(H3K4me3_intersect_infile, sep="\t", 
                                         names=["chr", "Start", "End", "Element.type", "ID", "peak_chr", "peak_start", "peak_end", "strand", "peak_ID", "depr"])
        H3K4me3_intersect = H3K4me3_intersect.drop(["strand", "depr"], axis=1)
        
        #element linked to H3K4me3 peaks they intersect with
        element_H3K4me3 = pd.merge(element, H3K4me3_intersect, on=['chr', 'Start', 'End', 'Element.type', 'ID'], how="left").fillna(0)
        #element_H3K4me3 = element_H3K4me3.drop_duplicates(subset=['chr', 'Start', 'End', 'Element.type', 'ID'], keep="first").fillna(0)
        element_H3K4me3[["peak_start", "peak_end", "peak_ID"]] = element_H3K4me3[["peak_start", "peak_end", "peak_ID"]].astype(int)
        
        # element connect to promoter through H3K4me3 peak
        element_promoter = pd.merge(element_H3K4me3, promoter_intersect_H3K4me3, on="peak_ID", how="left")
        element_promoter["H3K4me3"] = np.nan
        element_promoter.loc[element_promoter["peak_ID"] != 0, "H3K4me3"] = "H3K4me3"
        element_promoter["promoter"] = element_promoter["gene_name"]
        element_promoter[["H3K4me3", "promoter"]] = element_promoter[["H3K4me3", "promoter"]].astype(str)
        element_promoter = element_promoter.groupby(["chr", "Start", "End", "Element.type", "ID"], as_index=False).agg(",".join).reset_index().drop(["index"], axis=1)
               
        
        intron_intersect = pd.read_csv(intron_intersect_infile, sep="\t",
                                       names=["chr", "Start", "End", "Element.type", "ID", "intron_chr", "intron_start", "intron_end", "tx_id"])
        
        intron_intersect["transcript_id"] = intron_intersect["tx_id"].str.extract(r'(.+?)\.\d+')
        intron_intersect_update = pd.merge(intron_intersect, hg38_tx2gene, on=["transcript_id"], how="left")
        intron_intersect_update = intron_intersect_update.drop_duplicates(subset=['chr', 'Start', 'End', 'Element.type', 'ID', "gene_id"])
        # group by element
        intron_intersect_update = intron_intersect_update.groupby(["chr", "Start", "End", "Element.type", "ID"], as_index=False).agg(",".join)
        intron_intersect_update = intron_intersect_update[["chr", "Start", "End", "Element.type", "ID", "gene_name"]]
        
        #element annotated with intron
        element_intron = pd.merge(element, intron_intersect_update, on=["chr", "Start", "End", "Element.type", "ID"], how="left").fillna(0)
        element_intron["intron"] = float('NaN')
        element_intron.loc[element_intron["gene_name"] != 0, "intron"] = "intron"
        element_intron["intron_gene"] = element_intron["gene_name"]
        element_intron = element_intron.drop(["gene_name"], axis=1)
                
        
        exon_intersect = pd.read_csv(exon_intersect_infile, sep="\t",
                                     names=["chr", "Start", "End", "Element.type", "ID", "exon_chr", "exon_start", "exon_end", "tx_id", "exon_id"])
        exon_intersect["transcript_id"] = exon_intersect["tx_id"].str.extract(r'(.+?)\.\d+')
        exon_intersect_update = pd.merge(exon_intersect, hg38_tx2gene, on=["transcript_id"], how="left")
        exon_intersect_update = exon_intersect_update.drop_duplicates(subset=['chr', 'Start', 'End', 'Element.type', 'ID', "gene_id"])
        
        exon_intersect_update = exon_intersect_update.groupby(["chr", "Start", "End", "Element.type", "ID"], as_index=False).agg(",".join)
        exon_intersect_update = exon_intersect_update[["chr", "Start", "End", "Element.type", "ID", "gene_name"]]
        
        element_exon = pd.merge(element, exon_intersect_update, on=["chr", "Start", "End", "Element.type", "ID"], how="left").fillna(0)
        element_exon["exon"] = float('NaN')
        element_exon.loc[element_exon["gene_name"] != 0, "exon"] = "exon"
        element_exon["exon_gene"] = element_exon["gene_name"]
        element_exon = element_exon.drop(["gene_name"], axis=1)
        
        
        ann_list = [element_promoter, element_intron, element_exon]
        element_anno_complete = reduce(lambda left,right: pd.merge(left, right, on=["chr", "Start", "End", "Element.type", "ID"], how="outer"), ann_list)
        element_anno_complete[["intron", "exon"]] = element_anno_complete[["intron", "exon"]].astype(str)
        element_anno_complete.loc[element_anno_complete["intron_gene"] == 0, 'intron_gene'] = "nan"
        element_anno_complete.loc[element_anno_complete["exon_gene"] == 0, 'exon_gene'] = "nan"
        element_anno_complete["intergenic"] = "nan"
        element_anno_complete.loc[(element_anno_complete["H3K4me3"] == "nan") & (element_anno_complete["intron"] == "nan") & (element_anno_complete["exon"] == "nan"), "intergenic"] = "intergenic"
        
        element_anno_complete = element_anno_complete.replace(to_replace ="nan", value ="") 
        # element_anno_complete.to_csv(r'C:\Users\libin\UCSF\MMR\manual_annotation_complete_0718.csv', sep=",", index=False, header=True)
        
        
        ### group intron, exon to genic
        ### set annotation priority to : promoter, genic, intergenic
        element_anno_complete["annotation"] = np.nan
        element_anno_complete["cognate_gene"] = np.nan
        
        element_anno_complete.loc[element_anno_complete["promoter"] != "", "annotation"] = "promoter"
        element_anno_complete["cognate_gene"] = (element_anno_complete["promoter"]).where(element_anno_complete["annotation"].str.contains("promoter"))
        element_anno_complete.loc[(element_anno_complete["promoter"] == "") & ((element_anno_complete["exon"] != "") | (element_anno_complete["intron"] != "")), "annotation"] = "genic"
        element_anno_complete["cognate_gene"] = np.where((element_anno_complete["annotation"] == "genic"),(element_anno_complete["exon_gene"]+" "+element_anno_complete["intron_gene"]) ,element_anno_complete["cognate_gene"])
        element_anno_complete.loc[element_anno_complete["intergenic"] != "", "annotation"] = "intergenic"
        
        element_anno_complete_update = element_anno_complete[['chr','Start','End','Element.type','ID','annotation','cognate_gene']]
        element_anno_complete_update.to_csv(r'C:\Users\libin\UCSF\MMR\manual_annotation_complete_update_{}_{}.csv'.format(source, date), sep=",", index=False, header=True)

#%%

if chipseeker:
        promoter_ann = pd.read_csv(r'C:\Users\libin\UCSF\MMR\promoter.MMR.chipseek.output.bed', sep="\t")
        promoter_ann = promoter_ann[['seqnames', 'start', 'end', 'Element.type', 'ID', 'annotation', 'geneId']]
        promoter_ann["annotation"].loc[(promoter_ann["annotation"] == "Distal Intergenic")] = "Intergenic"
        promoter_ann["annotation_update"] = (promoter_ann["annotation"].astype(str) + " (" + promoter_ann["geneId"].astype(str) + ")").where(promoter_ann["annotation"].str.contains("Promoter"))
        promoter_ann["annotation_update"].fillna(promoter_ann["annotation"], inplace=True)
        promoter_ann = promoter_ann.drop(['annotation', 'geneId'], axis=1)
        promoter_ann["annotation"] = promoter_ann["annotation_update"].str.extract(r'([A-Za-z]+)(?=\s|$)')
        promoter_ann["gene_id"] = promoter_ann["annotation_update"].str.extract(r'\S+\s\(.*(ENSG.*?)\.')
        promoter_ann = pd.merge(promoter_ann, hg38_geneID2Name, on=["gene_id"] , how="left")
        promoter_ann = promoter_ann.drop(["annotation_update", "gene_type"], axis=1)
        
        #%%
        exon_ann = pd.read_csv(r'C:\Users\libin\UCSF\MMR\exon.MMR.chipseek.output.bed', sep="\t")
        exon_ann = exon_ann[['seqnames', 'start', 'end', 'Element.type', 'ID', 'annotation', 'geneId']]
        exon_ann["annotation"].loc[(exon_ann["annotation"] == "Distal Intergenic")] = "Intergenic"
        exon_ann["annotation_update"] = (exon_ann["annotation"].astype(str) + " (" + exon_ann["geneId"].astype(str) + ")").where(exon_ann["annotation"].str.contains("Promoter"))
        exon_ann["annotation_update"].fillna(exon_ann["annotation"], inplace=True)
        exon_ann = exon_ann.drop(['annotation', 'geneId'], axis=1)
        exon_ann["annotation"] = exon_ann["annotation_update"].str.extract(r'([A-Za-z]+)(?=\s|$)')
        exon_ann["gene_id"] = exon_ann["annotation_update"].str.extract(r'\S+\s\(.*(ENSG.*?)\.')
        exon_ann = pd.merge(exon_ann, hg38_geneID2Name, on=["gene_id"] , how="left")
        exon_ann = exon_ann.drop(["annotation_update", "gene_type"], axis=1)
        
        #%%
        intron_ann = pd.read_csv(r'C:\Users\libin\UCSF\MMR\intron.MMR.chipseek.output.bed', sep="\t")
        intron_ann = intron_ann[['seqnames', 'start', 'end', 'Element.type', 'ID', 'annotation', 'geneId']]
        intron_ann["annotation"].loc[(intron_ann["annotation"] == "Distal Intergenic")] = "Intergenic"
        intron_ann["annotation_update"] = (intron_ann["annotation"].astype(str) + " (" + intron_ann["geneId"].astype(str) + ")").where(intron_ann["annotation"].str.contains("Promoter"))
        intron_ann["annotation_update"].fillna(intron_ann["annotation"], inplace=True)
        intron_ann = intron_ann.drop(['annotation', 'geneId'], axis=1)
        intron_ann["annotation"] = intron_ann["annotation_update"].str.extract(r'([A-Za-z]+)(?=\s|$)')
        intron_ann["gene_id"] = intron_ann["annotation_update"].str.extract(r'\S+\s\(.*(ENSG.*?)\.')
        intron_ann = pd.merge(intron_ann, hg38_geneID2Name, on=["gene_id"] , how="left")
        intron_ann = intron_ann.drop(["annotation_update", "gene_type"], axis=1)
        
        #%%
        intergenic_ann = pd.read_csv(r'C:\Users\libin\UCSF\MMR\intergenic.MMR.chipseek.output.bed', sep="\t")
        intergenic_ann = intergenic_ann[['seqnames', 'start', 'end', 'Element.type', 'ID', 'annotation', 'geneId']]
        intergenic_ann["annotation"].loc[(intergenic_ann["annotation"] == "Distal Intergenic")] = "Intergenic"
        intergenic_ann["annotation_update"] = (intergenic_ann["annotation"].astype(str) + " (" + intergenic_ann["geneId"].astype(str) + ")").where(intergenic_ann["annotation"].str.contains("Promoter"))
        intergenic_ann["annotation_update"].fillna(intergenic_ann["annotation"], inplace=True)
        intergenic_ann = intergenic_ann.drop(['annotation', 'geneId'], axis=1)
        intergenic_ann["annotation"] = intergenic_ann["annotation_update"].str.extract(r'([A-Za-z]+)(?=\s|$)')
        intergenic_ann["gene_id"] = intergenic_ann["annotation_update"].str.extract(r'\S+\s\(.*(ENSG.*?)\.')
        intergenic_ann = pd.merge(intergenic_ann, hg38_geneID2Name, on=["gene_id"] , how="left")
        intergenic_ann = intergenic_ann.drop(["annotation_update", "gene_type"], axis=1)
        
        #%%
        ann_list = [promoter_ann, exon_ann, intron_ann, intergenic_ann]
        annotation_combined = reduce(lambda left,right: pd.merge(left, right, on=['seqnames', 'start', 'end', 'Element.type', 'ID', "annotation", "gene_id", "gene_name"], how="outer"), ann_list)
        annotation_combined[["gene_id", "gene_name"]] = annotation_combined[["gene_id", "gene_name"]].astype(str)
        
        #%%
        annotation_combined = annotation_combined.groupby(['seqnames', 'start', 'end', 'Element.type', 'ID']) \
                                                        .agg(",".join).reset_index()
        #%%
        annotation_split =  annotation_combined["annotation"].str.split(",").apply(pd.Series).add_prefix('annotation_')                                                
        gene_id_split =  annotation_combined["gene_id"].str.split(",").apply(pd.Series).add_prefix('gene_id_')                                                
        gene_name_split = annotation_combined["gene_name"].str.split(",").apply(pd.Series).add_prefix('gene_name_')
        
        #%%
        annotation_final = pd.concat([annotation_combined, annotation_split, gene_id_split, gene_name_split], axis=1)
        
        annotation_final = annotation_final[
        ['seqnames',
         'start',
         'end',
         'Element.type',
         'ID',
         'annotation_0',
         'gene_id_0',
         'gene_name_0',
         'annotation_1',
         'gene_id_1',
         'gene_name_1',
         'annotation_2',
         'gene_id_2',
         'gene_name_2']].fillna("")
        #%%
        annotation_final.to_csv(r'C:\Users\libin\UCSF\MMR\annotation_complete.csv', sep=",", index=False, header=True)