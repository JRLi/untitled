import pandas as pd
import os
import json
import numpy as np
import scipy.stats as st
ra_dir = 'E:/StringTemp/ra'


def open_df(in_path, direct='n'):
    fpath, fname = os.path.split(in_path)
    fbase, fext = os.path.splitext(fname)
    df = pd.read_csv(in_path, index_col=0) if fext == '.csv' else pd.read_table(in_path, index_col=0)
    if direct == 't':
        df = df.transpose()
    return df, fbase


def meta():
    df1, d1b = open_df(os.path.join(ra_dir, 'meta_s.txt'))
    df2, d2b = open_df(os.path.join(ra_dir, 'meta2_s.txt'))
    print(df1.shape, df2.shape)
    df3 = pd.concat([df1, df2], 1)
    df3 = df3.replace('#NULL!', np.nan)
    df3.to_csv(os.path.join(ra_dir, 'meta_test2.csv'), na_rep='NaN')


def manifest():
    with open(os.path.join(ra_dir, 'ra-assemb-manifest_181')) as in_181, \
            open(os.path.join(ra_dir, 'ra-assemb-manifest_all_reads')) as in_all, \
            open(os.path.join(ra_dir, 'meta_test2.txt')) as in_meta, open(os.path.join(ra_dir, 'man_3'), 'w') as out_f:
        next(in_181)
        out_f.write(next(in_all))
        next(in_meta)
        l_181 = [x.split(',')[0] for x in in_181]
        l_meta = [x.split('\t')[0] for x in in_meta]
        for line in in_all:
            lf = line.split(',')
            if lf[0] not in l_181 and lf[0] in l_meta:
                out_f.write(line)


def scipy_ttest_ind(s1, s2, var):
    return st.ttest_ind(s1, s2, equal_var=var)


def quantileNormalize(df_input):
    df = df_input.copy()
    #compute rank
    dic = {}
    for col in df:
        dic.update({col : sorted(df[col])})
    sorted_df = pd.DataFrame(dic)
    rank = sorted_df.mean(axis = 1).tolist()
    #sort
    for col in df:
        t = np.searchsorted(np.sort(df[col]), df[col])
        df[col] = [rank[i] for i in t]
    return df


def percentile(df_input, per_th = 10):
    df = df_input.copy()
    for i in range(0, len(df.columns)):
        qv = np.percentile(df.iloc[:, i], per_th)
        print(qv)
        df.iloc[:, i][df.iloc[:, i] < qv] = qv
        # print(df.iloc[:, i])
    return df


def df_mean_index(df_input):
    # gp1 = df1.groupby(df1.index)
    # df1 = gp1.mean()
    gp = df_input.groupby(df_input.index)
    return gp.mean()


def df_mean_columns(df_input):
    gp = df_input.groupby(df_input.columns, axis=1)
    return gp.mean()


l1 = [5.66, 5.78, 5.62, 5.65, 5.71, 5.73, 5.66, 5.71, 5.71, 5.62, 5.74, 5.68, 5.71, 5.83, 5.74, 5.70, 5.88, 5.79, 5.71, 5.72]
l2 = [5.58, 5.59, 5.66, 5.67, 5.67, 5.62, 5.63, 5.62, 5.55, 5.63, 5.58, 5.66, 5.68, 5.59, 5.65, 5.73, 5.71, 5.79, 5.67, 5.62]
na1 = np.array(l1)
na2 = np.array(l2)
na3 = np.log2(na1 + 0.01)
na4 = np.log2(na2 + 0.01)
print(scipy_ttest_ind(na1, na2, True))
print(scipy_ttest_ind(na1, na2, False))
print(scipy_ttest_ind(na3, na4, True))
print(scipy_ttest_ind(na3, na4, False))


dict1 = {'a': 3, 'b': 1, 'c': 2, 'd': 0}
dict1 = sorted(dict1, key=dict1.get, reverse=True)
print(dict1)

r_p = 'E:/StringTemp/Project_Rice/'
f_i = 'nm_df129_imputation.csv'
f_i2 = 'nm_df129.csv'
df, df_b = open_df(os.path.join(r_p, f_i))
df2, df_b2 = open_df(os.path.join(r_p, f_i2))
df.drop(['type (H)', 'waxy (H)'], axis=1, inplace=True)
df2.drop(['type (H)', 'waxy (H)'], axis=1, inplace=True)
df_t = df.iloc[:, 924:]
df_f = df.iloc[:, :924]
print(df_t.shape, df_f.shape)
df3 = df2.copy()
print(df_t.values.shape[0])
df2 = df2.dropna()
print(df2.shape)
df3 = df3.fillna(0)
print(df3.shape)
df3 = df3.dropna()
print(df3.shape)

with open('E:/StringTemp/lncRNA/pclnc_1rt0.5_2rt3.5_ct2_s1000.json', 'r') as json_file:
    subMatrix_list = json.load(json_file)
    i = 0
    for mod_dict in subMatrix_list:
        i += 1
        print('[mod{}]'.format(i))
        print('rows1', mod_dict.get('rows1'), sep='\n')
        print('rows2', mod_dict.get('rows2'), sep='\n')
        print('columns', mod_dict.get('columns'), sep='\n')

a = np.array([10])
b = np.array(10) -1
print(a)
print(b)
print(len(a))
print(type(b))
print(isinstance(10, int))

dfa, dfab = open_df(os.path.join(r_p, 'test.csv'))
print(dfa)
print(dfa.mean())
print(dfa.mean(1))
print(df_mean_index(dfa))