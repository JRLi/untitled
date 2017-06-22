#!/usr/bin/env python
import pandas as pd
import numpy as np
out_dir = 'LINCS_115209_split/'

df1 = pd.read_table('GSE70138_Level4_ZSVCINF_n115209x22268_20151231_annotation.txt', index_col=0)
df2 = pd.read_table('GSE70138_Level4_ZSVCINF_n115209x22268_20151231_all.txt', index_col=0)
print('df1 shape:', df1.shape)
print('df2 shape:', df2.shape)

len_count = 0
with open('cell_line_115209') as cl:
    for cell in cl:
        cell_line = cell.replace('\n', '')
        ii = np.where(df1.values == cell_line)[1]
        len_count += len(ii)
        print(cell_line, ii, len(ii))
        df3 = df2.iloc[:,ii]
        df3.to_csv('{}GSE70138_Level4_ZSVCINF_{}.csv'.format(out_dir, cell_line))

print('all len count:', len_count)