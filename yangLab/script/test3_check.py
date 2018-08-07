# -*- coding: utf-8 -*-
# The script is used to check the correlation between RPKM of two answers.
"""
Created on Thu Aug  2 16:09:39 2018

@author: liuzhen
"""

import numpy as np
import pandas as pd
from scipy import stats

result_path1 = '/Users/liuzhen/intern/data/test/RPKM_gene_dropscafd_result.txt'
result_path2 = '/Users/liuzhen/intern/data/test/zhn_RPKM.txt'

df1 = pd.read_table(result_path1, sep='\t', header=0)
df2 = pd.read_table(result_path2, sep='\t', header=None, names=['gene', 'chrom', 'RPKM'])
df2 = df2.sort_values(by=['gene'])
# df2.to_csv('/Users/liuzhen/intern/data/test/zhn_RPKM_sorted.txt', sep='\t', index=None)

df1.index = df1['gene']
df2.index = df2['gene']
df2_duplicates = df2[df2.duplicated('gene', keep=False)]
len_df2_dup = len(df2_duplicates.drop_duplicates('gene', keep='first'))
len_df2_all_gene = len(df2.drop_duplicates('gene', keep='first'))
print('Df2 duplicates gene number:', len_df2_dup)
print('Df2 duplicates gene percentage:', len_df2_dup/len_df2_all_gene)

# df2_duplicates.to_csv('/Users/liuzhen/intern/data/test/zhn_RPKM_dup.txt', sep='\t', index=None)

df2_rm_dup = pd.DataFrame(df2.drop_duplicates('gene', keep=False))
df1['RPKM2'] = df2_rm_dup['RPKM']
df1_uniq_gene = sum(df1.isna()['RPKM2'])
common_gene = sum(~df1.isna()['RPKM2'])
print('Df1 unique gene number:', df1_uniq_gene)
print('Common gene number:', common_gene)

df2_rm_dup['RPKM0'] = df1['RPKM']
df2_uniq_gene = sum(df2_rm_dup.isna()['RPKM0'])
common_gene = sum(~df2_rm_dup.isna()['RPKM0'])
print('Df2 remove duplicates unique gene number:', df2_uniq_gene)
print('Common gene number:', common_gene)

df1 = df1.dropna()
result_df = df1.loc[:, ['gene', 'RPKM', 'RPKM2']]
print(stats.pearsonr(result_df['RPKM'], result_df['RPKM2']))
# result_df.to_csv('/Users/liuzhen/intern/data/test/RPKM_check.txt', sep='\t', index=None)
