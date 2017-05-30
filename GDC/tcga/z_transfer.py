#!/usr/bin/env python
import sys
import os
import pandas as pd
import scipy.stats.mstats as mt


def openDF(in_path, direct = 'f'):
    fpath, fname = os.path.split(in_path)
    fbase, fext = os.path.splitext(fname)
    df = pd.read_csv(in_path, index_col=0) if fext == '.csv' else pd.read_table(in_path, index_col=0)
    if direct == 't':
        df = df.transpose()
    return df, fbase


def scipy_z_transfer(df_input):
    df = df_input.copy()
    for k in range(len(df.columns)):
        df.iloc[:, k] = mt.zscore(df.iloc[:, k], ddof=1)
    return df


def df_column_stat(df_input):
    df = df_input.copy()
    return df.mean(), df.std(), df.min(), df.max()


def main():
    if len(sys.argv) != 2:
        sys.exit(0)
    df1, df_base = openDF(sys.argv[1])
    print(df1.shape)
    df1 = scipy_z_transfer(df1)
    print(df1.shape)
    avg, std, minV, maxV = df_column_stat(df1)
    # df4 = pd.DataFrame(columns=['avg', 'std', 'max', 'min'])
    df2 = pd.DataFrame()
    df2['avg'] = avg
    df2['std'] = std
    df2['max'] = maxV
    df2['min'] = minV
    df1.to_csv('{}_zScore.csv'.format(df_base))
    df2.to_csv('{}_summary.csv'.format(df_base))

if __name__ == '__main__':
    main()