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


def df_stat(df_input, direct='col_wise'):
    df = df_input.copy()
    df2 = pd.DataFrame()
    axis_dir = 0 if direct == 'col_wise' else 1
    df2['avg'] = df.mean(axis_dir)
    df2['std'] = df.std(axis_dir)
    df2['max'] = df.max(axis_dir)
    df2['min'] = df.min(axis_dir)
    return df2


def main(argv=None):
    if argv is None:
        argv = args_parse()
    print(argv)
    df1, df_base = openDF(argv.inputDF)
    print(df1.shape)
    df1 = scipy_z_transfer(df1, argv.direct)
    trans = '_col_wiseZ' if argv.direct == 'c' else '_row_wiseZ'
    print(df1.shape)
    df_cwise = df_stat(df1, 'col_wise')
    df_iwise = df_stat(df1, 'index_wise')
    df1.to_csv('{}{}.csv'.format(df_base, trans))
    df_cwise.to_csv('Summary_col_{}{}.csv'.format(df_base, trans))
    df_iwise.to_csv('Summary_row_{}{}.csv'.format(df_base, trans))

if __name__ == '__main__':
    sys.exit(main())