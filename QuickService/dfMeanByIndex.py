#!/usr/bin/env python3.6
import pandas as pd
import sys
import os


def openDF(in_path, direct='f'):
    fpath, fname = os.path.split(in_path)
    fbase, fext = os.path.splitext(fname)
    df = pd.read_csv(in_path, index_col=0) if fext == '.csv' else pd.read_table(in_path, index_col=0)
    if direct == 't':
        df = df.transpose()
    return df, fbase


def df_mean_index(df_input):
    gp = df_input.groupby(df_input.index)
    print('mean group by index.')
    return gp.mean()


def main():
    df1, fbase = openDF(sys.argv[1])
    print(fbase, 'default shape:', df1.shape)
    df1 = df_mean_index(df1)
    print(fbase, 'mean shape:', df1.shape)
    df1.to_csv('{}_m.csv'.format(fbase))

if __name__ == '__main__':
    main()