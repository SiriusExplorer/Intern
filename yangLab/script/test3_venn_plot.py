# The script is used to draw a venn plot for the test3 RPKM result between mine and the answer.

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

result_path1 = r'..\data\result\test3\RPKM_gene_dpsfd_avglen_result.txt'
result_path2 = r'..\data\result\test3\zhn_RPKM.txt'

df1 = pd.read_table(result_path1, sep='\t', header=0)
df2 = pd.read_table(result_path2, sep='\t', header=None, names=['gene', 'chrom', 'RPKM'])
df2.drop_duplicates('gene', keep='first')
df2 = df2.sort_values(by=['gene'])
df2.index = range(len(df2))

plt.figure()
venn2([set(df1['gene']), set(df2['gene'])], ('new', 'old'))
plt.show()



