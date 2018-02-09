import pandas as pd
import os
from collections import defaultdict


def open_df(in_path, direct='n'):
    fpath, fname = os.path.split(in_path)
    fbase, fext = os.path.splitext(fname)
    df = pd.read_csv(in_path, index_col=0) if fext == '.csv' else pd.read_table(in_path, index_col=0)
    print('{}: Transpose: {}'.format(fbase, direct))
    if direct == 't':
        df = df.transpose()
    return df, fbase



def main():
    df_label, label_base = open_df('D:/Project/circ_miRNA/out_annotation/lab_TCGA-CHOL.csv', 't')
    print(df_label)
    print(df_label.iloc[0:, 1])

if __name__ == '__main__':
    main()