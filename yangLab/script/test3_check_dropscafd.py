# -*- coding: utf-8 -*-
# The script is used to check the correlation between RPKM of two answers.
"""
Created on Thu Aug  2 16:09:39 2018

@author: liuzhen
"""

import numpy as np
import pandas as pd
from scipy import stats

result_path1 = '/Users/liuzhen/intern/data/test/result/test3/RPKM_gene_dpsfd_avglen_result.txt'
result_path2 = '/Users/liuzhen/intern/data/test/zhn_RPKM.txt'

df1 = pd.read_table(result_path1, sep='\t', header=0)
df2 = pd.read_table(result_path2, sep='\t', header=None, names=['gene', 'chrom', 'RPKM'])
df2 = df2.sort_values(by=['gene'])
df2.index = range(len(df2))
# df2.to_csv('/Users/liuzhen/intern/data/test/zhn_RPKM_sorted.txt', sep='\t', index=None)

bool_list = []
bool_list2 = []
for i in range(len(df2)):
    if '_' in df2.loc[i, 'chrom']:
        bool_list.append(False)
        bool_list2.append(True)
    else:
        bool_list.append(True)
        bool_list2.append(False)
df2_dropscafd = df2[bool_list]
df2_dropgene = df2[bool_list2]
print('The gene number of dropped scaffold:', len(df2_dropgene.drop_duplicates('gene', keep='first')))

df1.index = df1['gene']
df2_dropscafd.index = df2_dropscafd['gene']
df2_dpsfd_duplicates = df2_dropscafd[df2_dropscafd.duplicated('gene', keep=False)]
len_df2_all_gene = len(df2.drop_duplicates('gene', keep='first'))
len_df2_dpsfd_dup = len(df2_dpsfd_duplicates.drop_duplicates('gene', keep='first'))
len_df2_dpsfd_all_gene = len(df2_dropscafd.drop_duplicates('gene', keep='first'))
print('The gene been totally dropped:', len_df2_all_gene-len_df2_dpsfd_all_gene)
print('Df2 drop scaffold duplicates gene number:', len_df2_dpsfd_dup)
print('Df2 drop scaffold duplicates gene percentage:', len_df2_dpsfd_dup/len_df2_dpsfd_all_gene)

# df2_dpsfd_duplicates.to_csv('/Users/liuzhen/intern/data/test/zhn_RPKM_dpsfd_dup.txt', sep='\t', index=None)
# df2_duplicates.to_csv('/Users/liuzhen/intern/data/test/zhn_RPKM_dup.txt', sep='\t', index=None)

df2_dpsfd_rm_dup = pd.DataFrame(df2_dropscafd.drop_duplicates('gene', keep='first'))
df1['RPKM2'] = df2_dpsfd_rm_dup['RPKM']
df1_uniq_gene = sum(df1.isna()['RPKM2'])
common_gene = sum(~df1.isna()['RPKM2'])
print('Df1 unique gene number:', df1_uniq_gene)
print('Common gene number:', common_gene)

df2_dpsfd_rm_dup['RPKM0'] = df1['RPKM']
df2_dropscafd = pd.DataFrame(df2_dropscafd)
df2_dropscafd['RPKM0'] = df1['RPKM']
df2_dropscafd = df2_dropscafd.dropna()
df2_dpsfd_dup = df2_dropscafd[df2_dropscafd.duplicated('gene', keep=False)]
df2_dpsfd_dup_gene = len(df2_dpsfd_dup.drop_duplicates('gene', keep='first'))
print('The number of duplicated genes in the common genes:', df2_dpsfd_dup_gene)
df2_uniq_gene = sum(df2_dpsfd_rm_dup.isna()['RPKM0'])
common_gene = sum(~df2_dpsfd_rm_dup.isna()['RPKM0'])
print('Df2 drop scaffold unique gene number:', df2_uniq_gene)
print('Common gene number:', common_gene)

df1 = df1.dropna()
result_df = df1.loc[:, ['gene', 'RPKM', 'RPKM2']]
print('PCC by genes:', stats.pearsonr(result_df['RPKM'], result_df['RPKM2']))
# result_df.to_csv('/Users/liuzhen/intern/data/test/RPKM_check.txt', sep='\t', index=None)

# df2_dropscafd = pd.DataFrame(df2_dropscafd)
# df2_dropscafd['RPKM0'] = df1['RPKM']
# df2_dropscafd = df2_dropscafd.dropna()
# result_df = df2_dropscafd.loc[:, ['gene', 'RPKM0', 'RPKM']]
# print(stats.pearsonr(result_df['RPKM0'], result_df['RPKM']))