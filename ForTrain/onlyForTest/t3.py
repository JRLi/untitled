import pandas as pd
import numpy as np
import datetime
import os
import dask.dataframe as dd
import scipy.stats as st
import scipy.stats.mstats as mt
from sklearn.preprocessing import StandardScaler, scale


def openDF(in_path, direct = 'f'):
    fpath, fname = os.path.split(in_path)
    fbase, fext = os.path.splitext(fname)
    df = pd.read_csv(in_path, index_col=0) if fext == '.csv' else pd.read_table(in_path, index_col=0)
    if direct == 't':
        df = df.transpose()
    return df, fbase


def scipy_z_transfer2(df_input, dof=1, direct='c'):
    df = df_input.copy()
    if direct == 'c':
        for k in range(len(df.columns)):
            df.iloc[:, k] = mt.zscore(df.iloc[:, k], ddof=dof)
    else:
        for k in range(len(df.index)):
            df.iloc[k, :] = mt.zscore(df.iloc[k, :], ddof=dof)
    return df


def df_stat(df_input, direct='col_wise'):
    df = df_input.copy()
    df2 = pd.DataFrame()
    axis_dir = 0 if direct == 'col_wise' else 1
    df2['avg'] = df.mean(axis_dir)
    df2['std'] = df.std(axis_dir)
    df2['max'] = df.max(axis_dir)
    df2['min'] = df.min(axis_dir)
    return df2


def s_top(series_input, top_n):
    ss = series_input.copy()
    ss.sort_values(inplace=True)
    ss = ss.iloc[range(-top_n, top_n)]
    return ss


def scipy_corr(s1, s2, corr_mode):
    ixs = s1.index.intersection(s2.index)
    if corr_mode == 'pearson':
        return st.pearsonr(s1[ixs], s2[ixs])
    elif corr_mode == 'kendall':
        return st.kendalltau(s1[ixs], s2[ixs])
    elif corr_mode == 'spearman':
        return st.spearmanr(s1[ixs], s2[ixs])


def corr_by_col_of_df(df_c, df_d, top, corr_m):
    dfc = pd.DataFrame(columns=df_c.columns, index=df_d.columns)
    dfp = pd.DataFrame(columns=df_c.columns, index=df_d.columns)
    print('target_shape:', dfc.shape)
    count_c, count_all = 0, 0
    for i in range(len(df_c.columns)):
        c1 = df_c.columns[i]
        s1 = df_c[c1]
        count_c += 1
        c_list, p_list = [], []
        for j in range(len(df_d.columns)):
            c2 = df_d.columns[j]
            count_all += 1
            s2 = df_d[c2] if top in (0, None) else s_top(df_d[c2], top)
            cor_e, p_value = scipy_corr(s1, s2, corr_m)
            c_list.append(cor_e)
            p_list.append(p_value)
        dfc[c1] = pd.Series(c_list).values
        dfp[c1] = np.array(p_list)
    print('Cells count: {}, Drugs count: {}'.format(count_c, count_all / count_c))
    return dfc, dfp

t1 = datetime.datetime.now()

df1, df_base = openDF('D://Project/drs/forTest9x7.txt')
ddf = dd.read_table('D://Project/RA/in/*')
print(df1)
print(ddf.head())

df2 = scipy_z_transfer2(df1)
df3 = scipy_z_transfer2(df1, 0)
df4 = df_stat(df2)
df5 = df_stat(df3)
print(df4)
print(df5)

sc = StandardScaler()
sc.fit(df1)
df6 = pd.DataFrame(sc.transform(df1), index=df1.index, columns=df1.columns)
df7 = df_stat(df6)
print(df7)
df8 = pd.DataFrame(scale(df1), index=df1.index, columns=df1.columns)
df9 = df_stat(df8)
print(df9)

df10 = pd.read_table('D://Project/drs/gdsc/LINCS/GSE70138_Level4_ZSVCINF_n115209x22268_20151231_annotation.txt', index_col=0)
print(df10.shape)
df10.T.to_csv('D://Project/drs/gdsc/LINCS/GSE70138_20151231_annotation.txt', sep='\t')
td = datetime.datetime.now() - t1
print("\t[Info] Spending time={0}!".format(td))