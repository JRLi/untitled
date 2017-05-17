#!/usr/bin/env python
import pandas as pd
import os
import sys


def openDF(in_path, direct = 'f'):
    fpath, fname = os.path.split(in_path)
    fbase, fext = os.path.splitext(fname)
    df = pd.read_csv(in_path, index_col=0) if fext == '.csv' else pd.read_table(in_path, index_col=0)
    if direct == 't':
        df = df.transpose()
    return df, fbase


def df_column_stat(df_input):
    df = df_input.copy()
    return df.mean(), df.std(), df.min(), df.max()


def main():
    if len(sys.argv) != 2:
        print('no input file')
        sys.exit(0)
    df_path = sys.argv[1]
    df1, df_base = openDF(df_path)  # if need transpose, set direct  to 't'
    avg, std, minV, maxV = df_column_stat(df1)
    df2 = pd.DataFrame(columns=['avg', 'std', 'max', 'min'])
    df2['avg'] = avg
    df2['std'] = std
    df2['max'] = maxV
    df2['min'] = minV
    df2.to_csv(df_base + '_columns_status.csv')

if __name__ == '__main__':
    main()