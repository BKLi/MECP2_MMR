# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 22:25:54 2020

@author: bingkun
"""

import pandas as pd

gene_list = pd.read_csv(r'C:\Users\libin\UCSF\MMR\refseq.gene.list.byTX', sep="\t", names=["chr", "start", "end", "name"])
gene_list["length"] = gene_list["end"] - gene_list["start"]
gene_list = gene_list.sort_values(by=['length'], ascending=False)
#gene_list = gene_list.replace({"NC_000002.12": "chr2", "NC_000003.12": "chr3", "NC_000007.14":"chr7", "NC_000020.11":"chr20", "NC_000023.11":"chrX"})
element_weak = pd.read_csv(r'C:\Users\libin\UCSF\MMR\enrichment\element_weak', sep="\t", names=["chr", "start", "end"])
element_strong = pd.read_csv(r'C:\Users\libin\UCSF\MMR\enrichment\element_strong', sep="\t", names=["chr", "start", "end"])
element_pos = pd.concat([element_weak, element_strong])

#%%
# HPRT1 chrX:134460164 - 134500668
chrX = gene_list[gene_list["chr"] == "chrX"]
HPRT1 = chrX[chrX["name"] == "HPRT1"].reset_index()
target_gene_start = HPRT1.loc[0, "start"]
target_gene_end = HPRT1.loc[0, "end"]

element_pos_chrX = element_pos[element_pos["chr"] == "chrX"]
gene_list_row = []
for index, row in element_pos_chrX.iterrows(): 
        ele_start = row["start"]
        ele_end = row["end"]
        if ele_end < target_gene_start:
                inter_gene_list = chrX[(chrX["start"] > ele_end) & (chrX["end"] < target_gene_start)]
        elif ele_start > target_gene_end:
                inter_gene_list = chrX[(chrX["start"] > target_gene_end) & (chrX["end"] < ele_start)]
        inter_gene_list = inter_gene_list.drop_duplicates(subset=["name"])
        print (inter_gene_list['name'].tolist())
        gene_list_row.append(inter_gene_list['name'].tolist())
        gene_list_row_num = [len(i) for i in gene_list_row]
element_pos_chrX.loc[:, "inter_genes"] = gene_list_row
element_pos_chrX.loc[:, "inter_genes_number"] = gene_list_row_num


#%%
# MLH1 chr3:35,993,000-38,051,000 
chr3 = gene_list[gene_list["chr"] == "chr3"]
MLH1 = chr3[chr3["name"] == "MLH1"].reset_index()
target_gene_start = MLH1.loc[0, "start"]
target_gene_end = MLH1.loc[0, "end"]

element_pos_chr3 = element_pos[element_pos["chr"] == "chr3"]
gene_list_row = []
for index, row in element_pos_chr3.iterrows(): 
        ele_start = row["start"]
        ele_end = row["end"]
        if ele_end < target_gene_start:
                inter_gene_list = chr3[(chr3["start"] > ele_end) & (chr3["end"] < target_gene_start)]
        elif ele_start > target_gene_end:
                inter_gene_list = chr3[(chr3["start"] > target_gene_end) & (chr3["end"] < ele_start)]
        inter_gene_list = inter_gene_list.drop_duplicates(subset=["name"])
        print (inter_gene_list['name'].tolist())
        gene_list_row.append(inter_gene_list['name'].tolist())
        gene_list_row_num = [len(i) for i in gene_list_row]
element_pos_chr3.loc[:, "inter_genes"] = gene_list_row
element_pos_chr3.loc[:, "inter_genes_number"] = gene_list_row_num


#%%
# PMS2 chr7:4,973,000-7,010,000 
chr7 = gene_list[gene_list["chr"] == "chr7"]
PMS2 = chr7[chr7["name"] == "PMS2"].reset_index()
target_gene_start = PMS2.loc[0, "start"]
target_gene_end = PMS2.loc[0, "end"]

element_pos_chr7 = element_pos[element_pos["chr"] == "chr7"]
gene_list_row = []
for index, row in element_pos_chr7.iterrows(): 
        ele_start = row["start"]
        ele_end = row["end"]
        if ele_end < target_gene_start:
                inter_gene_list = chr7[(chr7["start"] > ele_end) & (chr7["end"] < target_gene_start)]
        elif ele_start > target_gene_end:
                inter_gene_list = chr7[(chr7["start"] > target_gene_end) & (chr7["end"] < ele_start)]
        inter_gene_list = inter_gene_list.drop_duplicates(subset=["name"])
        print (inter_gene_list['name'].tolist())
        gene_list_row.append(inter_gene_list['name'].tolist())
        gene_list_row_num = [len(i) for i in gene_list_row]
element_pos_chr7.loc[:, "inter_genes"] = gene_list_row
element_pos_chr7.loc[:, "inter_genes_number"] = gene_list_row_num

#%%
# PCNA chr20:4,114,000-6,127,000
chr20 = gene_list[gene_list["chr"] == "chr20"]
PCNA = chr20[chr20["name"] == "PCNA"].reset_index()
target_gene_start = PCNA.loc[0, "start"]
target_gene_end = PCNA.loc[0, "end"]

element_pos_chr20 = element_pos[element_pos["chr"] == "chr20"]
gene_list_row = []
for index, row in element_pos_chr20.iterrows(): 
        ele_start = row["start"]
        ele_end = row["end"]
        if ele_end < target_gene_start:
                inter_gene_list = chr20[(chr20["start"] > ele_end) & (chr20["end"] < target_gene_start)]
        elif ele_start > target_gene_end:
                inter_gene_list = chr20[(chr20["start"] > target_gene_end) & (chr20["end"] < ele_start)]
        inter_gene_list = inter_gene_list.drop_duplicates(subset=["name"])
        print (inter_gene_list['name'].tolist())
        gene_list_row.append(inter_gene_list['name'].tolist())
        gene_list_row_num = [len(i) for i in gene_list_row]
element_pos_chr20.loc[:, "inter_genes"] = gene_list_row
element_pos_chr20.loc[:, "inter_genes_number"] = gene_list_row_num

#%%
inter_gene_number = pd.concat([element_pos_chr3, element_pos_chr7, element_pos_chr20, element_pos_chrX])
inter_gene_number.to_csv(r'C:\Users\libin\UCSF\MMR\number_between_gene_and_enhancer.csv', sep=",", index=False)

#%%    
# MSH2-MSH6 chr2:46,400,000-48,810,000 
#chr2 = gene_list[gene_list["chr"] == "chr2"]
#MSH2_MSH6_list = chr2[(chr2["start"] > 46400000) & (chr2["end"] < 48810000 )]
#MSH2_MSH6_list["length"] = MSH2_MSH6_list["end"] - MSH2_MSH6_list["start"]
#MSH2_MSH6_list.to_csv(r'C:\Users\libin\UCSF\MMR\MSH2_MSH6.target.geneList', sep="\t", index=False)
