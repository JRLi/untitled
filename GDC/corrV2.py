#!/usr/bin/env python
import sys
import os
import numpy as np
import pandas as pd
import argparse
import datetime
import scipy.stats as st

use_message = '''
    Need Python3 and numpy, pandas, scipy; or Anaconda.
    Usage: corrV3.py [-t, -c] df1 df2
    Example: corrV3.py -t 250 -c s df1.csv df2.txt 
    OR: corrV3.py df1 df2, will use all genes and pearson correlation.
'''


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


def args_parse():
    parser = argparse.ArgumentParser(description=use_message)
    parser.add_argument('-t', '--top', type=int, help='Top and Down -t number genes(only select the second file); '
                                                      'if set 0 equal No use -t')
    parser.add_argument('-c', '--corr', type=str, choices=['p', 'k', 's'],
                        default='p', help='Correlation method, p: pearson, k: kendall, s: spearman; default is p')
    parser.add_argument('-d', '--direct', choices=['n', 't'],
                        default=['n', 'n'], nargs=2, help="n is normal, t is transpose; default is n n")
    parser.add_argument("-m", "--mean", action="store_true", help="row-wise mean for the same index")
    parser.add_argument('pairs', nargs=2, help="cell line (1) and drug (2) expression profile")
    args = parser.parse_args()
    return args


def prepare_output_dir(output_dir):
    if os.path.exists(output_dir):
        pass
    else:
        os.mkdir(output_dir)


def openDF(in_path, direct):
    fpath, fname = os.path.split(in_path)
    fbase, fext = os.path.splitext(fname)
    df = pd.read_csv(in_path, index_col=0) if fext == '.csv' else pd.read_table(in_path, index_col=0)
    print('{}: Transpose: {}'.format(fbase, direct))
    if direct == 't':
        df = df.transpose()
    return df, fbase


def df_mean_index(df_input):
    gp = df_input.groupby(df_input.index)
    return gp.mean()


def s_top(series_input, top_n):
    ss = series_input.copy()
    ss.sort_values(inplace=True)
    ss = ss.iloc[range(-top_n, top_n)]
    return ss


def scipy_corr(s1, s2, corr_mode):
    ixs = s1.index.intersection(s2.index)
    if corr_mode == 'pearson':
        return st.pearsonr(s1[ixs], s2[ixs])
    elif corr_mode == 'kendall':
        return st.kendalltau(s1[ixs], s2[ixs])
    elif corr_mode == 'spearman':
        return st.spearmanr(s1[ixs], s2[ixs])


def corr_by_col_of_df(df_cell, df_drug, top, corr_m):
    dfc = pd.DataFrame(columns=df_drug.columns, index=df_cell.columns)
    dfp = pd.DataFrame(columns=df_drug.columns, index=df_cell.columns)
    print('target_shape:', dfc.shape)

    count_drug, count_all = 0, 0
    for i in range(len(df_drug.columns)):
        col_drug = df_drug.columns[i]
        series_drug = df_drug[col_drug] if top in (0, None) else s_top(df_drug[col_drug], top)
        count_drug += 1
        corr_list, pvalue_list = [], []
        for j in range(len(df_cell.columns)):
            col_cell = df_cell.columns[j]
            series_cell = df_cell[col_cell]
            count_all += 1
            cor_e, p_value = scipy_corr(series_cell, series_drug, corr_m)
            corr_list.append(cor_e)
            pvalue_list.append(p_value)
        dfc[col_drug] = pd.Series(corr_list).values
        dfp[col_drug] = np.array(pvalue_list)

    print('Drugs count: {}, Cells count: {}'.format(count_drug, count_all / count_drug))
    return dfc, dfp


def main(argv=None):
    try:
        if argv is None:
            time_1 = datetime.datetime.now()
            print('start time:', str(time_1))
            argv = args_parse()
            df_cells, cells_base = openDF(argv.pairs[0], argv.direct[0])
            df_drugs, drugs_base = openDF(argv.pairs[1], argv.direct[1])
            print('file1: ' + cells_base, 'file2: ' + drugs_base, sep='\n')
            print('df1 rows: {}\tdf1 columns: {}'.format(df_cells.shape[0], df_cells.shape[1]))
            print('df2 rows: {}\tdf2 columns: {}'.format(df_drugs.shape[0], df_drugs.shape[1]))

            prepare_output_dir('./correlation_output')
            prepare_output_dir('./p_value_output')

            if argv.mean:
                df_cells, df_drugs = df_mean_index(df_cells), df_mean_index(df_drugs)
                print('df_cells.shape after mean normalization:', df_cells.shape)
                print('df_drugs.shape after mean normalization:', df_drugs.shape)

            df_cells = df_cells.loc[(df_cells!=0).any(1), (df_cells!=0).any(0)]
            df_drugs = df_drugs.loc[(df_drugs!=0).any(1), (df_drugs!=0).any(0)]
            print('df_cells.shape after remove all zero rows and columns:', df_cells.shape)
            print('df_drugs.shape after remove all zero rows and columns:', df_drugs.shape)
            ixs = df_cells.index.intersection(df_drugs.index)
            print('The numbers of intersection index:', len(ixs))

            method = {'p': 'pearson', 'k': 'kendall', 's': 'spearman'}.get(argv.corr)
            print('correlation mode:', method)
            top_suffix = 'all' if argv.top in (0, None) else 'top' + str(argv.top)
            print('Use Top:', top_suffix)
            df5c, df5p = corr_by_col_of_df(df_cells, df_drugs, argv.top, method)
            time_2 = datetime.datetime.now()
            print('after correlation:', str(time_2))

            df5c.to_csv('./correlation_output/Corr_{}_{}_{}_{}.csv'.format(cells_base, drugs_base, method, top_suffix))
            df5p.to_csv('./p_value_output/P_value_{}_{}_{}_{}.csv'.format(cells_base, drugs_base, method, top_suffix))
            time_3 = datetime.datetime.now()
            print('Finished time:', str(time_3))
            print('All used time:', str(time_3 - time_1))
        print('Done.\n')

    except Usage as err:
        print(sys.stderr, err.msg)
        print(sys.stderr, "Terminated, for help use -h or --help")
        return 2

if __name__ == "__main__":
    sys.exit(main())