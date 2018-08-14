# The script is used to draw a boxplot for the PCC result from test1

import numpy as np
import pandas as pd
from scipy import stats
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib import axes

def draw_heatmap(data, xlabels, ylabels):
    figure, ax = plt.subplots(facecolor='w')
    ax.set_yticks(range(len(ylabels)))
    ax.set_yticklabels(ylabels)
    ax.set_xticks(range(len(xlabels)))
    ax.set_xticklabels(xlabels)
    map = ax.imshow(data, interpolation='nearest', cmap=cm.Blues, vmin=data.min(), vmax=data.max())
    plt.colorbar(mappable=map, cax=None, ax=None, shrink=0.5)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    plt.show()
    return 0

data_df = pd.read_table(r'..\data\raw_data\test1\pcc.txt', sep='\t', header=0)
sample_list = list(data_df.columns[1:])
sample_list.sort()
pcc_matrix = []
for sample1 in sample_list:
    pcc_list = []
    for sample2 in sample_list:
        pcc_list.append(stats.pearsonr(data_df[sample1], data_df[sample2])[0])
    pcc_matrix.append(pcc_list)
pcc_npm = np.array(pcc_matrix)

draw_heatmap(pcc_npm, sample_list, sample_list)