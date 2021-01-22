# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 13:52:28 2020

@author: bingkun
!!! This one does not work -- 10.19.2020
"""
import pandas as pd
import numpy as np
import itertools

binsize = 5000
chroms = ["chr21"]
MACS2_PATH = r'C:\Users\libin\UCSF\MMR\plac\combined.1e-3_merged.narrowPeak'
MACS2_full = pd.read_csv(MACS2_PATH, sep='\t', skip_blank_lines=True,comment='#',header=None, usecols=[0,1,2])
MACS2_full.columns = ['chr','start','end']
MACS2_full['start_bin'] = np.floor(MACS2_full['start']/binsize)*binsize
MACS2_full['end_bin'] = np.ceil(MACS2_full['end']/binsize)*binsize
for CHR in chroms:
        print('doing chromosome ',CHR,'\n')
        #handling MACS2 peaks
        print('-- handling MACS2 peaks')
        MACS2 = MACS2_full[MACS2_full['chr'] == CHR].copy()
        MACS2['start_bin'] = np.floor(MACS2['start']/binsize)
        MACS2['end_bin'] = np.ceil(MACS2['end']/binsize)
        #perform this hack becasue apply returns wrong data type in some rare case
        specialCase = False
        if MACS2.iloc[0]['end_bin'] - MACS2.iloc[0]['start_bin'] == MACS2.shape[1] - 1:
            MACS2.iloc[0,MACS2.columns.get_loc('start_bin')] = MACS2.iloc[0]['start_bin'] - 1
            specialCase = True
        MACS2_peak_ranges = MACS2.apply(lambda row: range(int(row['start_bin']),int(row['end_bin'])), axis=1).values.tolist()
        MACS2_peak_ranges_list = set(itertools.chain.from_iterable(MACS2_peak_ranges))
        if specialCase:
            MACS2_peak_ranges_list.remove(MACS2.iloc[0]['start_bin'])
