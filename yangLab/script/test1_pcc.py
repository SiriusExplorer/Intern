# The script is used to calculate the pearson correlation coefficient of test1 data

import numpy as np
import pandas as pd
from scipy import stats

def cal_PCC(list1, list2):
    n1 = np.array(list1)
    n2 = np.array(list2)
    n = len(n1)
    n1a = np.average(n1)
    n2a = np.average(n2)
    PCC = (np.sum(n1*n2)-n*n1a*n2a) / (np.sqrt(np.sum(n1*n1)-n*n1a**2)*np.sqrt(np.sum(n2*n2)-n*n2a**2))
    return PCC

data_df = pd.read_table('/home/liuzhen2018/data/test/pcc.txt', sep='\t', header=0)
sample_list = list(data_df.columns[1:])
result_df = pd.DataFrame()
i = 0
for sample1 in sample_list:
    for sample2 in sample_list[sample_list.index(sample1)+1:]:
        result_df.loc[i, 'sample1'] = sample1
        result_df.loc[i, 'sample2'] = sample2
        result_df.loc[i, 'PCC'] = cal_PCC(data_df[sample1], data_df[sample2])
        # result_df.loc[i, 'PCC'] = stats.pearsonr(data_df[sample1], data_df[sample2])[0]
        i += 1
result_df.to_csv('/home/liuzhen2018/data/test/pcc_result.txt', sep='\t', index=False, header=False)


