# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 16:38:58 2020

@author: bingkun
@project: MMR

plot ranks of candidate genes (in quantile)
*input: TMM-normalized RPKM
*output: scatterplot
"""
#%%
import pandas as pd
import seaborn as sns
from scipy import stats
import matplotlib
import matplotlib.pyplot as plt 
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
pd.set_option('display.float_format', lambda x: '%.5f' % x)
pd.options.display.precision = 5

#%%
sample_name = ["LM124", "LM134"]
gene_list = ["HPRT1", "MSH2", "MSH6", "MLH1", "PMS2",
             "PCNA", "SOX2", "NANOG", "POU5F1", "KLF4"]

#%%
# read in TMM-normalized RPKM expression values
expression = pd.read_csv(r'C:\Users\libin\UCSF\MECP2\RNASeq\TMM\TMM_RPKM_DG', sep="\t")
expression = expression.reset_index()
expression["gene_id"] = expression["index"].str.extract(r'(.+?)\.\d+')
expression = expression.drop(["index"], axis=1)
# 
expression[sample_name] = expression[sample_name].fillna(0)
expression["average_expression"] = expression[sample_name].mean(axis=1)
#%%
#
hg38_geneID2geneName = pd.read_csv(r'C:\Users\libin\UCSF\hg38_general\gencode.v32.geneID2geneNameType', sep='\t')

# calc percentile rank
expression_merged = pd.merge(expression, hg38_geneID2geneName, how="inner").reset_index()
expression_merged = expression_merged.drop(["index"], axis=1)
# https://stackoverflow.com/questions/44211653/calculate-percentile-for-every-value-in-a-column-of-dataframe
expression_merged["rank"] = expression_merged["average_expression"].rank(method='min', pct=True).round(4)#.astype("int")
expression_merged["average_expression"] = expression_merged["average_expression"].round(3)
sns.distplot(expression_merged["average_expression"], hist=False, kde=True)
expression_merged["average_expression"].median()
expression_merged["average_expression"].mean()

expression_selected = expression_merged[expression_merged["gene_name"].isin(gene_list)]

expression_merged.to_csv(r'C:\Users\libin\UCSF\MECP2\RNASeq\TMM\average_expression.csv', sep="," , index=False, header=True)
#%%
# make quantile plot
fig = plt.figure(figsize=(10,10))
sns.scatterplot(x=expression_selected["rank"], y=expression_selected["average_expression"])


#%% ---- test code ----
#expression_merged["pct_rank_{}".format(sample_name[0])] = expression_merged["{}".format(sample_name[0])].rank(method='min', pct=True).round(3)#.astype("int")
#expression_merged["pct_rank_{}".format(sample_name[1])] = expression_merged["{}".format(sample_name[1])].rank(method="min", pct=True).round(3)#.astype('int')
#expression_merged["pct_rank_{}".format(sample_name[0])] = [stats.percentileofscore(expression_merged["{}".format(sample_name[0])].values, i, kind='strict') for i in expression_merged["{}".format(sample_name[0])].values]
#expression_merged["pct_rank_{}".format(sample_name[1])] = [stats.percentileofscore(expression_merged["{}".format(sample_name[1])].values, i, kind='strict') for i in expression_merged["{}".format(sample_name[1])].values]

#expression_merged["pct_rank_{}".format(sample_name[0])] = expression_merged["pct_rank_{}".format(sample_name[0])]/expression_merged.shape[0]
#expression_merged["pct_rank_{}".format(sample_name[1])] = expression_merged["pct_rank_{}".format(sample_name[1])]/expression_merged.shape[0]

#expression_merged[sample_name] = expression_merged[sample_name].astype("int64")

#expression_selected.to_csv(r'C:\Users\libin\UCSF\MECP2\RNASeq\TMM\test', sep="\t", index=False)
