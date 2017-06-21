import pandas as pd
import numpy as np
import os
import scipy.stats.mstats as mt
import json

sample = 'abcde'
lf = sample.split('|')
print(lf)
s1 = sample[:sample.rfind('d')]
print(s1)
s2 = sample[:sample.rfind('z')]
print(sample.rfind('z'))
print(s2)


def openDF(in_path, direct = 'f'):
    fpath, fname = os.path.split(in_path)
    fbase, fext = os.path.splitext(fname)
    df = pd.read_csv(in_path, index_col=0) if fext == '.csv' else pd.read_table(in_path, index_col=0)
    if direct == 't':
        df = df.transpose()
    return df, fbase


def z_transfer(df_input, mode='n'):
    df = df_input.copy()
    print('z score transfer.')
    # mean and std track column from row index 0 to end, so axis use default (index 0).
    avg = df.mean()
    std = df.std()
    for k in range(len(df.columns)):
        df.iloc[:, k] = (df.iloc[:, k] - avg[k]) / std[k]
    if mode == 'n':
        return df
    else:
        return -df


def scipy_z_transfer(df_input):
    df = df_input.copy()
    for k in range(len(df.columns)):
        df.iloc[:, k] = mt.zscore(df.iloc[:, k], ddof=1)
    return df


def scipy_z_transfer2(df_input, direct='c'):
    df = df_input.copy()
    if direct == 'c':
        for k in range(len(df.columns)):
            df.iloc[:, k] = mt.zscore(df.iloc[:, k], ddof=1)
    else:
        for k in range(len(df.index)):
            df.iloc[k, :] = mt.zscore(df.iloc[k, :], ddof=1)
    return df


def df_column_stat(df_input):
    df = df_input.copy()
    return df.mean(), df.std(), df.min(), df.max()


def df_index_stat(df_input):
    df = df_input.copy()
    return df.mean(1), df.std(1), df.min(1), df.max(1)

df1, df_base = openDF('D://Project/drs/forTest9x7.txt')
print(df1)
df2 = z_transfer(df1)
df3 = scipy_z_transfer(df1)
df6 = scipy_z_transfer2(df1, 'i')
print(df2)
print(df3)
print(df6)

avg, std, minV, maxV = df_column_stat(df2)
avgi, stdi, minVi, maxVi = df_index_stat(df6)
#df4 = pd.DataFrame(columns=['avg', 'std', 'max', 'min'])
df4 = pd.DataFrame()
df4['avg'] = avg
df4['std'] = std
df4['max'] = maxV
df4['min'] = minV
print(df4)
df5 = pd.DataFrame()
df5['avg'] = avgi
df5['std'] = stdi
df5['max'] = maxVi
df5['min'] = minVi
print(df5)

with open('D://Project/drs/GDC/projects.json') as json_file, open('D://Project/drs/GDC/projects.txt', 'w') as outFile:
    dict_list = json.load(json_file)
    for mod_dict in dict_list:
        pid = mod_dict.get('project_id')
        psite = ''.join(mod_dict.get('primary_site'))
        dtype = ''.join(mod_dict.get('disease_type'))
        print('{}\t{}\t{}'.format(pid, dtype, psite))
        outFile.write('{}\t{}\t{}\n'.format(pid, dtype, psite))