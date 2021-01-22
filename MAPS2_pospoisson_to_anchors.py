# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 16:01:17 2020

@author: bingkun
* get the full list of anchors from MAPS output
"""

import pandas as pd
import numpy as np

chrom_list = ["chr3", "chr7", "chr20"]
for chrom in chrom_list:
        #chrom = "chr21"
        AND_interactions = pd.read_csv('reg_raw.{}.combined.5k.and.MAPS2_pospoisson'.format(chrom), delim_whitespace=True)
        # anchor list part one
        AND_anchor1 = AND_interactions[["chr", "bin1_mid", "X1D_peak_bin1"]]
        AND_anchor2 = AND_interactions[["chr", "bin2_mid", "X1D_peak_bin2"]]
        
        XOR_interactions = pd.read_csv('reg_raw.{}.combined.5k.xor.MAPS2_pospoisson'.format(chrom), delim_whitespace=True)
        # anchor list part two
        xor_bin1_anchor = XOR_interactions[(XOR_interactions["X1D_peak_bin1"] == 1) & (XOR_interactions["X1D_peak_bin2"] == 0)]        
        xor_bin1_anchor = xor_bin1_anchor[["chr", "bin1_mid", "X1D_peak_bin1"]]
        # anchor list part three
        xor_bin2_anchor = XOR_interactions[(XOR_interactions["X1D_peak_bin1"] == 0) & (XOR_interactions["X1D_peak_bin2"] == 1)]       
        xor_bin2_anchor = xor_bin2_anchor[["chr", "bin2_mid", "X1D_peak_bin2"]]

        anchor_list = pd.DataFrame(np.concatenate([AND_anchor1, AND_anchor2, xor_bin1_anchor, xor_bin2_anchor]), columns=AND_anchor1.columns)
        anchor_list = anchor_list.drop_duplicates()
        
        anchor_list["start"] = anchor_list["bin1_mid"].astype('int64')
        #anchor_list["start"] = anchor_list["start"].astype('int64')
        anchor_list["end"] = anchor_list["bin1_mid"] + 5000
        anchor_list["end"] = anchor_list["end"].astype('int64')
        anchor_list = anchor_list[["chr", "start", "end"]]
        anchor_list["anchor_ID"] = ["{}_anchor_{}".format(chrom,i) for i in range(anchor_list.shape[0])]
        anchor_list.to_csv('{}_anchor_list'.format(chrom),sep="\t", index=False, header=False)
        
        print ('number of anchors on {} in total: {}'.format(chrom, anchor_list.shape[0]))
