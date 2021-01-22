# -*- coding: utf-8 -*-
"""
Created on Sun Oct 18 23:55:38 2020

@author: libin
"""

import pandas as pd
#%%
chrom_list = ["chr3", "chr7", "chr20"]
for chrom in chrom_list:
        matrix_file = '{}.normalized.hic.input'.format(chrom)
        matrix = pd.read_csv(matrix_file, sep="\t", names=["chr", "start1", "end1", "start2", "end2", "count"])
        # remove self-interacting bins
        matrix = matrix[matrix["start1"] != matrix["start2"]]
        #matrix = matrix[matrix["count"] != "NaN"]
        matrix["ID"] = ["{}_inter_{}".format(chrom, i) for i in range(matrix.shape[0])]
        
        matrix_bin1 = pd.DataFrame()
        matrix_bin1["chr"] = matrix["chr"]
        matrix_bin1["start"] = matrix["start1"]
        matrix_bin1["end"] = matrix["end1"]
        matrix_bin1["ID"] = matrix["ID"]
        matrix_bin1["count"] = matrix["count"]        
        matrix_bin1.to_csv(r'{}.normalized.hic.bin1'.format(chrom), sep="\t", header=False, index=False)
        
        matrix_bin2 = pd.DataFrame()
        matrix_bin2["chr"] = matrix["chr"]
        matrix_bin2["start"] = matrix["start2"]
        matrix_bin2["end"] = matrix["end2"]
        matrix_bin2["ID"] = matrix["ID"]
        matrix_bin2["count"] = matrix["count"]        
        matrix_bin2.to_csv(r'{}.normalized.hic.bin2'.format(chrom), sep="\t", header=False, index=False)
        
        matrix_inter = matrix[["chr", "start1", "start2", "ID", "count"]]
        matrix_inter.to_csv(r'{}.inter'.format(chrom), sep="\t", header=True, index=False)
