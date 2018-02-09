#!/usr/bin/env python3.6
import pandas as pd
import os
from sklearn.model_selection import LeaveOneOut
from sklearn import svm
f_dir = 'D:/Project/lncRNA/exoRbase'
n_p = 'N'
d_p = 'HCC'


def open_df(in_path, direct='n'):
    fpath, fname = os.path.split(in_path)
    fbase, fext = os.path.splitext(fname)
    df = pd.read_csv(in_path, index_col=0) if fext == '.csv' else pd.read_table(in_path, index_col=0)
    if direct == 't':
        df = df.transpose()
    return df, fbase


def main():
    df1, dfb = open_df(os.path.join(f_dir, 'mRNA_TPM_10genes.csv'))
    n_c = [col for col in df1 if col.startswith(n_p) or col.startswith(d_p)]
    df1 = df1[n_c]
    print(df1)
    for i in range(len(df1)):
        run = df1[0:1]
        print(run)
        df1 = df1[1:]


if __name__ == '__main__':
    main()
