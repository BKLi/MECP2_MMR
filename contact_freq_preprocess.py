# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 01:24:43 2020

@author: bingkun
"""

import pandas as pd
#%%
chrom_list = ["chrX", "chr3", "chr7", "chr20"]
res = "5000"
for chrom in chrom_list:
        matrix_file = r'C:\Users\libin\UCSF\MMR\hic\{}.{}.txt'.format(chrom, res)
        matrix = pd.read_csv(matrix_file, sep="\t", names=["start1", "start2", "count"])
        # remove self-interacting bins
        matrix = matrix[matrix["start1"] != matrix["start2"]]
        #matrix = matrix[matrix["count"] != "NaN"]
        matrix["chr"] = chrom
        matrix["ID"] = ["{}_inter_{}".format(chrom, i) for i in range(matrix.shape[0])]
        
        matrix_bin1 = pd.DataFrame()
        matrix_bin1["chr"] = matrix["chr"]
        matrix_bin1["start"] = matrix["start1"]
        matrix_bin1["end"] = matrix_bin1["start"] + int(res)
        matrix_bin1["ID"] = matrix["ID"]
        matrix_bin1["count"] = matrix["count"]        
        matrix_bin1.to_csv(r'C:\Users\libin\UCSF\MMR\hic\processed\{}_{}.bin1'.format(chrom, res), sep="\t", header=False, index=False)
        
        matrix_bin2 = pd.DataFrame()
        matrix_bin2["chr"] = matrix["chr"]
        matrix_bin2["start"] = matrix["start2"]
        matrix_bin2["end"] = matrix_bin2["start"] + int(res)
        matrix_bin2["ID"] = matrix["ID"]
        matrix_bin2["count"] = matrix["count"]        
        matrix_bin2.to_csv(r'C:\Users\libin\UCSF\MMR\hic\processed\{}_{}.bin2'.format(chrom, res), sep="\t", header=False, index=False)

        matrix_inter = matrix[["chr", "start1", "start2", "ID", "count"]]
        matrix_inter.to_csv(r'C:\Users\libin\UCSF\MMR\hic\processed\{}_{}.inter'.format(chrom, res), sep="\t", header=True, index=False)
#%%
element_strong = pd.read_csv(r'C:\Users\libin\UCSF\MMR\enrichment\element_strong', sep="\t", names=["chr", "start", "end"])
element_weak = pd.read_csv(r'C:\Users\libin\UCSF\MMR\enrichment\element_weak', sep="\t", names=["chr", "start", "end"])
element_negative = pd.read_csv(r'C:\Users\libin\UCSF\MMR\enrichment\element_negative', sep="\t", names=["chr", "start", "end"])

promoter = pd.read_csv(r'C:\Users\libin\UCSF\hg38_general\refseq\refseq.hg38.byTX.1kb.promoter.correct', sep="\t", names=["chr", "start", "end", "tx_id", "gene_name"])
gene_list = ["MLH1", "PMS2", "PCNA", "HPRT1"]
promoter_select = promoter[promoter["gene_name"].isin(gene_list)]

element_strong["ID"] = ["strong_{}".format(i) for i in range(element_strong.shape[0])]
element_weak["ID"] = ["weak_{}".format(i) for i in range(element_weak.shape[0])]

element_positive = pd.concat([element_strong, element_weak])

element_promoter = pd.merge(element_positive, promoter_select, on=["chr"], how="inner")
element_promoter["dist"] = abs(element_promoter["start_y"] - element_promoter["start_x"])

for gene in gene_list:
        promoter_in_list = promoter[promoter["gene_name"] == gene]
        chrom_gene = promoter_in_list["chr"].tolist()[0]
        promoter_in_list.to_csv(r'C:\Users\libin\UCSF\MMR\hic\analysis\promoter_{}'.format(chrom_gene), sep="\t", header=False, index=False)
        
        element_strong_in_list = element_strong[element_strong["chr"] == chrom_gene]
        element_strong_in_list.to_csv(r'C:\Users\libin\UCSF\MMR\hic\analysis\element_strong_{}'.format(chrom_gene), sep="\t", header=False, index=False)
        element_weak_in_list = element_weak[element_weak["chr"] == chrom_gene]
        element_weak_in_list.to_csv(r'C:\Users\libin\UCSF\MMR\hic\analysis\element_weak_{}'.format(chrom_gene), sep="\t", header=False, index=False)
        element_negative_in_list = element_negative[element_negative["chr"] == chrom_gene]
        element_negative_in_list.to_csv(r'C:\Users\libin\UCSF\MMR\hic\analysis\element_negative_{}'.format(chrom_gene), sep="\t", header=False, index=False)

#%%
chrX_inter = pd.read_csv(r'C:\Users\libin\UCSF\MMR\hic\processed\chrX_5000.inter'.format(chrom, res), sep="\t")
chrX_inter_select = chrX_inter[chrX_inter["start1"] == 135095000]
chrX_inter_select_2 = chrX_inter[chrX_inter["start2"] == 135095000]
chrX_inter_select_3 = chrX_inter[chrX_inter["start1"] == 134455000]

chr7_inter = pd.read_csv(r'C:\Users\libin\UCSF\MMR\hic\processed\chr7_5000.inter'.format(chrom, res), sep="\t")
chr7_inter_select = chr7_inter[chr7_inter["start1"] == 5410000]
chr7_inter_select_2 = chr7_inter[chr7_inter["start2"] == 5410000]
chr7_inter_select_3 = chr7_inter[chr7_inter["start2"] == 5990000]
