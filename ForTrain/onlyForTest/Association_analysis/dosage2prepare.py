#!/usr/bin/env python
import os
import pandas as pd
f_dos = '/mount/ictr1/chenglab/jli/Projects/RNA_expression/GTEx/dosage-out1-rm-iloc'
f_tpm = '/mount/ictr1/chenglab/jli/Projects/RNA_expression/GTEx/gtex_v6p_RPKM_Lung_donor_mean_gn'


def df_open(path_in):
    p_d, f_n = os.path.split(path_in)
    f_p, f_s = os.path.splitext(f_n)
    df = pd.read_csv(path_in, index_col=0) if f_s == '.csv' else pd.read_table(path_in, index_col=0)
    return df, p_d, f_p


def main():
    df1, dir_p1, file_p1 = df_open(f_dos)
    df1 = df1.iloc[:, 5:]
    df1.columns = ['-'.join(x.split('-')[:2]) for x in df1.columns]
    df1.sort_index(axis=1, inplace=True)
    # df1.to_csv(os.path.join(dir_p, '{}-df'.format(file_p)), sep='\t')

    df2, dir_p2, file_p2 = df_open(f_tpm)
    df2.sort_index(axis=1, inplace=True)
    idx = df2.columns.intersection(df1.columns)
    df1 = df1[idx]
    df2 = df2[idx]
    print(df1.shape, df2.shape)
    df1.to_csv(os.path.join(dir_p1, '{}-inter-df'.format(file_p1)), sep='\t')
    df2.to_csv(os.path.join(dir_p2, '{}_inter'.format(file_p2)), sep='\t')


if __name__ == '__main__':
    main()
