#!/usr/bin/env python
import pandas as pd
import os

fileCheck = 'snv_indel'
stringCheck = 'HASH'


def openDF(in_path, direct='f'):
    fpath, fname = os.path.split(in_path)
    fbase, fext = os.path.splitext(fname)
    df = pd.read_csv(in_path, index_col=0) if fext == '.csv' else pd.read_table(in_path, index_col=0)
    if direct == 't':
        df = df.transpose()
    return df, fbase


def main():
    linc_list = os.listdir('.')
    for lincs in linc_list:
        if lincs.startswith(fileCheck):
            df1, df1_base = openDF(lincs)
            print(df1_base, df1.shape)
            cols = [c for c in df1.columns if stringCheck not in c]
            df1 = df1[cols]
            print('after remove target:', df1.shape)
            df1.to_csv(df1_base + '.csv')

if __name__ == '__main__':
    main()