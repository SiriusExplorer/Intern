# The script is used to draw a boxplot for the PCC results from test1

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def print_data(mod, df):
    print("{0} PCC median:".format(mod), np.median(df['PCC']))
    print("{0} PCC 1/4:".format(mod), np.percentile(df['PCC'], 25))
    print("{0} PCC 3/4:".format(mod), np.percentile(df['PCC'], 75))
    return 0

column_list = ['sample1', 'sample2', 'PCC']
data_df = pd.read_table(r'..\data\raw_data\test1\pcc.txt', sep='\t', header=0)
sample_list = list(data_df.columns[1:])
hand_df = pd.read_table(r'..\data\result\test1\pcc_result.txt', sep='\t', header=None, names=column_list)
stats_df = pd.read_table(r'..\data\result\test1\pcc_result_stats.txt', sep='\t', header=None, names=column_list)

plt.figure(figsize=(6, 3))
ax = data_df.boxplot(column=sample_list, grid=False)
plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
plt.title('Sample Value Range')
plt.show()

# plt.figure()
# plt.subplot(131)
#
# integrate_df = pd.DataFrame(hand_df)
# integrate_df['PCC_scipy'] = stats_df['PCC']
# integrate_df.boxplot(column=['PCC', 'PCC_scipy'], grid=False, widths=0.5)
# plt.ylim(0.8, 1)
# plt.title('Sample Correlation')
# plt.show()
# print_data('hand', hand_df)
# print_data('scipy', stats_df)