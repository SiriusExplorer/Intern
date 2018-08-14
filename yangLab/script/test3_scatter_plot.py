# The script is used to draw a scatter plot for the test3 RPKM result between mine and the answer.

import pandas as pd
import numpy as np
from scipy import stats
import statsmodels.api as sm
import matplotlib.pyplot as plt

def print_data(mod, df):
    print("{0} PCC 1/4:".format(mod), np.percentile(df['RPKM'], 25))
    print("{0} PCC median:".format(mod), np.median(df['RPKM']))
    print("{0} PCC 3/4:".format(mod), np.percentile(df['RPKM'], 75))
    return 0

result_path1 = r'..\data\result\test3\RPKM_gene_dpsfd_avglen_result.txt'
result_path2 = r'..\data\result\test3\zhn_RPKM.txt'

df1 = pd.read_table(result_path1, sep='\t', header=0)
df2 = pd.read_table(result_path2, sep='\t', header=None, names=['gene', 'chrom', 'RPKM'])

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

df1.index = df1['gene']
df2_dropscafd.index = df2_dropscafd['gene']
df2_dpsfd_rm_dup = pd.DataFrame(df2_dropscafd.drop_duplicates('gene', keep='first'))
df1['RPKM2'] = df2_dpsfd_rm_dup['RPKM']
df2_dropscafd = pd.DataFrame(df2_dropscafd)
df2_dropscafd['RPKM0'] = df1['RPKM']
df1 = df1.dropna()
df2_dropscafd = df2_dropscafd.dropna()
df1.index = range(len(df1))
df2_dropscafd.index = range(len(df2_dropscafd))

print(df1[df1['RPKM'] == max(df1['RPKM'])])

fig, ax = plt.subplots()
ver_n1 = np.array(df1['RPKM'])
ver_n2 = np.array(df2_dropscafd['RPKM'])
ver_n1.shape = (-1, 1)
ver_n2.shape = (-1, 1)
ax.boxplot([ver_n1, ver_n2], showfliers=False)
ax.set_xticklabels(['new', 'origin'])
ax.set_title('RPKM Destribution')
plt.show()

print_data('new', df1)
print_data('origin', df2_dropscafd)

# integrate_df = df1.loc[:, ['gene', 'RPKM', 'RPKM2']]
# integrate_df.to_csv(r'..\data\result\test3\check_rpkm.txt', sep='\t', index=None)

# y = np.array(df1['RPKM2'])
# x = np.array(df1['RPKM'])
# X = sm.add_constant(x)
# lm = sm.OLS(y, X).fit()
# print(lm.summary())
# fig, ax = plt.subplots()
# ax.scatter(df1['RPKM'], df1['RPKM2'], s=5, label='data')
# ax.plot(x, lm.fittedvalues, 'r--', linewidth=1, label='OLS')
# ax.legend(loc='best')
# ax.set_title('Linear Regression between RPKMs (remove duplicated genes)')
# ax.set_xlabel('new RPKM')
# ax.set_ylabel('original RPKM')
# plt.show()

