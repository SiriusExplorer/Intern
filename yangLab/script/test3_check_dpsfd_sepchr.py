# -*- coding: utf-8 -*-
# The script is used to check the correlation between RPKM of two answers.
# The RPKM is calculated based on genes at each chromosome
"""
Created on Thu Aug  2 16:09:39 2018

@author: liuzhen
"""

import numpy as np
import pandas as pd
from scipy import stats

result_path1 = '/Users/liuzhen/intern/data/test/result/test3/RPKM_gene_dpsfd_avglen_sepchr_result.txt'
result_path2 = '/Users/liuzhen/intern/data/test/zhn_RPKM.txt'

df1 = pd.read_table(result_path1, sep='\t', header=0)
df2 = pd.read_table(result_path2, sep='\t', header=None, names=['gene', 'chrom', 'RPKM'])
df2 = df2.sort_values(by=['gene'])
df2.index = range(len(df2))
# df2.to_csv('/Users/liuzhen/intern/data/test/zhn_RPKM_sorted.txt', sep='\t', index=None)

df1.index = [df1['gene'], df1['chrom']]
df2.index = [df2['gene'], df2['chrom']]

df1['RPKM2'] = df2['RPKM']
df1 = df1.dropna()
df1.index = range(len(df1))


result_df = df1.loc[:, ['gene', 'RPKM', 'RPKM2']]
result_df_dup = result_df[result_df.duplicated('gene', keep=False)]
print('The number of duplicated genes repeats:', len(result_df_dup))
result_df_dup_gene = len(result_df_dup.drop_duplicates('gene', keep='first'))
print('The number of duplicated genes in the common genes:', result_df_dup_gene)
result_df.drop_duplicates('gene', keep=False, inplace=True)
print(stats.pearsonr(result_df['RPKM'], result_df['RPKM2']))
# result_df.to_csv('/Users/liuzhen/intern/data/test/RPKM_check.txt', sep='\t', index=None)
