#!/usr/bin/env python
import pandas as pd
import numpy as np
import os
''''a = {}
b = set()
print(type(a), type(b))

a = {1, 3, 4}
b = {1, 4 ,5}
print(type(a))
print(a.union(b))
d = a.intersection(b)
print(a)
a.update(b)
print(a)

a = [1, 2, 3]
a.append('5')
print(a)
print(b)    # return None to b
print('xxxxx')
a = 'abcde'
b = a[::-1]
a = a.replace('a', 'e').replace('e', 'c').replace('c', 't')
print(b)
print(a)
ab = []

aa = ['1','2','3','50','4']
for i in range(1,len(aa)):
    if int(aa[i]) - int(aa[i-1]) ==1:
        ab.append(aa[i])
print(ab)
'''''
seq_string = 'ACTGatcgnN'
seq_r_string = seq_string[::-1]
complement_dict = {'a':'t','t':'a','c':'g','g':'c','A':'T','T':'A','C':'G','G':'C'}
print('\t', '[INFO]',seq_r_string)
seq_rc_list =[]
for a in seq_r_string:
    seq_rc_list.append(complement_dict.get(a))
    if a == 'A':
        seq_rc_list.append('T')
    elif a == 'C':
        seq_rc_list.append('G')
    elif a == 'T':
        seq_rc_list.append('A')
    elif a == 'G':
        seq_rc_list.append('C')
    elif a == 'a':
        seq_rc_list.append('t')
    else:
        seq_rc_list.append(a)
print(seq_rc_list)
#seqRC = "".join(seq_rc)
seq_rc_string = ''
for a in seq_rc_list:
    seq_rc_string += '1'
print(seq_rc_string)

st1 = '445\t333\tfff'
lf = st1.split('\t')
a, b, c = lf
print(a, b, c)
'''
cell_line = 'JURKAT'
df1 = pd.read_table('GSE70138_Level4_ZSVCINF_n115209x22268_20151231_annotation.txt', index_col=0)
print('df1_shape:', df1.shape)
print('ii')
ii = np.where(df1.values == cell_line)[1]
print(ii)
jj = np.where(df1.values == cell_line)
print('jj')
print(jj)
'''
def openDF(in_path, direct = 'f'):
    fpath, fname = os.path.split(in_path)
    fbase, fext = os.path.splitext(fname)
    df = pd.read_csv(in_path, index_col=0) if fext == '.csv' else pd.read_table(in_path, index_col=0)
    if direct == 't':
        df = df.transpose()
    return df, fbase

df1, df_base = openDF('D://Project/drs/forTest9x7.txt')
avg = df1.mean()
std = df1.std()
minV = df1.min()
maxV = df1.max()
print(avg)
print(std)
print(minV)
print(maxV)

df2 = pd.DataFrame(columns=['avg', 'std', 'max', 'min'])
df2['avg'] = avg
df2['std'] = std
df2['max'] = maxV
df2['min'] = minV
print(df2)

df1 = pd.read_table('D://Project/drs/GDC/bb.txt', index_col=0, header=None, names=['bb'])
print(df1)
df2 = pd.read_table('D://Project/drs/GDC/aa.txt', index_col=0, header=None, names=['aa'])
df3 = pd.DataFrame()
df3['bb'] = df1['bb']
df3['aa'] = df2['aa']
#df3.rename_axis('aaaa')
print(df3)
print(df3.shape)
print(df3.index)