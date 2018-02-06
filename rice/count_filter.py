#!/usr/bin/env python3.6
import os
import numpy as np
import pandas as pd
from sklearn.preprocessing import Imputer


def df_open(path_in, direct='n'):
    d_p, d_n = os.path.split(path_in)
    f_p, f_t = os.path.splitext(d_n)
    df1 = pd.read_csv(path_in, index_col=0) if f_t == '.csv' else pd.read_table(path_in, index_col=0)
    if direct != 'n':
        df1 = df1.T
    return df1, d_p, f_p


def sk_imputation(df_in):
    df = df_in.copy()
    imr = Imputer(missing_values='NaN', strategy='mean', axis=0)
    imr = imr.fit(df)
    arr = imr.transform(df)
    return arr


def rpm(df_input):
    df = df_input.copy()
    return (df / df.sum()) * 1000000


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


def df_sum(df_in, direct='n'):
    ax = 0 if direct == 'n' else 1
    return df_in.sum(ax)


def d_filter(df_in, sum_ss, number):
    return df_in.loc[:, sum_ss >= number]


def phenotype_concat(df_in, f_name):
    df_p, p_path, p_pre = df_open('E:/StringTemp/Project_Rice/phenotype.csv', 't')
    df_p = df_p.add_suffix('p')
    df_x = pd.concat([df_in, df_p], join='inner')
    df_x.T.to_csv(os.path.join(p_path, '{}.csv'.format(f_name)))
    arr_i = sk_imputation(df_x)
    df_i = pd.DataFrame(arr_i, index=df_x.index, columns=df_x.columns)
    df_i.T.to_csv(os.path.join(p_path, '{}_i.csv'.format(f_name)))


def main():
    df_w, w_path, w_pre = df_open('E:/StringTemp/Project_Rice/rice_nov_mir_wm.csv')
    df_n, n_path, n_pre = df_open('E:/StringTemp/Project_Rice/rice_nov_mir_nm.csv')
    print(df_n.shape)
    for i in range(20000, 100001, 10000):
        print('K'.join(str(i).rsplit('000', 1)))
        df_wp = d_filter(df_w, df_sum(df_w), i)
        df_np = d_filter(df_n, df_sum(df_n), i)
        print(df_wp.shape)
        print(df_np.shape)
        if i == 50000:
            df_w50k_rpm = rpm(df_wp)
            df_n50k_rpm = rpm(df_np)
            df_w50k_qn = quantileNormalize(df_wp)
            df_n50k_qn = quantileNormalize(df_np)
            phenotype_concat(df_w50k_rpm, 'wm_df50k_r')
            phenotype_concat(df_n50k_rpm, 'nm_df50k_r')
            phenotype_concat(df_w50k_qn, 'wm_df50k_q')
            phenotype_concat(df_n50k_qn, 'nm_df50k_q')


if __name__ == '__main__':
    main()