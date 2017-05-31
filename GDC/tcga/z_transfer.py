#!/usr/bin/env python
import sys
import os
import argparse
import pandas as pd
import scipy.stats.mstats as mt


def args_parse():
    parser = argparse.ArgumentParser(description='z_transfer.py')
    parser.add_argument('-d', '--direct', choices=['c', 'i'],
                        default='c', help="standardizing within a column (c) or a index (i); default is c")
    parser.add_argument('-i', '--inputDF', help='data frame before split')
    args = parser.parse_args()
    return args


def openDF(in_path, direct = 'f'):
    fpath, fname = os.path.split(in_path)
    fbase, fext = os.path.splitext(fname)
    df = pd.read_csv(in_path, index_col=0) if fext == '.csv' else pd.read_table(in_path, index_col=0)
    if direct == 't':
        df = df.transpose()
    return df, fbase


def scipy_z_transfer(df_input, direct):
    df = df_input.copy()
    if direct == 'c':
        for k in range(len(df.columns)):
            df.iloc[:, k] = mt.zscore(df.iloc[:, k], ddof=1)
    else:
        for k in range(len(df.index)):
            df.iloc[k, :] = mt.zscore(df.iloc[k, :], ddof=1)
    return df


def df_column_wise_stat(df_input):
    df = df_input.copy()
    return df.mean(), df.std(), df.max(), df.min()


def df_index_wise_stat(df_input):
    df = df_input.copy()
    return df.mean(1), df.std(1),  df.max(1), df.min(1)


def create_stat_df(avg_series, std_series, max_series, min_series):
    df = pd.DataFrame()
    df['avg'] = avg_series
    df['std'] = std_series
    df['max'] = max_series
    df['min'] = min_series
    return df


def main(argv=None):
    if argv is None:
        argv = args_parse()
    print(argv)
    df1, df_base = openDF(argv.inputDF)
    print(df1.shape)
    df1 = scipy_z_transfer(df1, argv.direct)
    trans = '_col_wiseZ' if argv.direct == 'c' else '_row_wiseZ'
    print(df1.shape)
    avg_c, std_c, max_c, min_c = df_column_wise_stat(df1)
    avg_i, std_i, max_i, min_i = df_index_wise_stat(df1)
    df_c = create_stat_df(avg_c, std_c, max_c, min_c)
    df_i = create_stat_df(avg_i, std_i, max_i, min_i)
    df1.to_csv('{}{}.csv'.format(df_base, trans))
    df_c.to_csv('{}_column_wise.csv'.format(df_base))
    df_i.to_csv('{}_index_wise.csv'.format(df_base))

if __name__ == '__main__':
    sys.exit(main())