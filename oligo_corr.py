#!/usr/bin/env python
# coding: utf-8

# In[11]:


import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


# In[2]:


raw_counts = pd.read_csv(r'C:\Users\libin\UCSF\MECP2\MMR_sub_readsCRISPRn_22264682.txt', sep='\t')


# In[5]:


corr_matrix = raw_counts.drop(["oligo"], axis=1).corr()


# In[13]:


plt.figure(figsize=(10,10))
sns.clustermap(corr_matrix)
plt.savefig(r'C:\Users\libin\UCSF\MECP2\MMR_sub_readsCRISPRn_22264682_corr.pdf', transparent=True)




