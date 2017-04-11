#!/usr/bin/env python
import sys
import os
import numpy as np
import pandas as pd
import argparse

use_message = '''
    Need Python3 and numpy, pandas, scipy.
    Usage: correlationOfDF.py df1 df2
'''


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


def args_parse():
    parser = argparse.ArgumentParser(description=use_message)
    parser.add_argument('-t', '--top', type=int, help='base pairs per line; if set 0 equal No use top')
    parser.add_argument('pairs', nargs=2, help="cell line (1) and drug (2) expression profile")
    args = parser.parse_args()
    return args


def openDF(in_path):
    fpath, fname = os.path.split(in_path)
    fbase, fext = os.path.splitext(fname)
    df = pd.read_table(in_path, index_col=0)
    return df, fbase


def df_mean_index(df_input):
    gp = df_input.groupby(df_input.index)
    return gp.mean()


def s_top(series_input, top_n):
    ss = series_input.copy()
    ss.sort_values(inplace=True)
    ss = ss.iloc[range(-top_n, top_n)]
    return ss


def main(argv=None):
    try:
        if argv is None:
            argv = args_parse()
            df_cells, cells_base = openDF(argv.pairs[0])
            df_drugs, drugs_base = openDF(argv.pairs[1])
            df_cells, df_drugs = df_mean_index(df_cells), df_mean_index(df_drugs)
            print('df_cells.shape:', df_cells.shape)
            print('df_drugs.shape:', df_drugs.shape)
            # df_cd = df_cells.merge(df_drugs,left_index=True, right_index=True, how='left').dropna()
            # create output data frames for corr and p-value
            df5c = pd.DataFrame(columns=df_cells.columns, index=df_drugs.columns)
            df5p = pd.DataFrame(columns=df_cells.columns, index=df_drugs.columns)
            print('target_shape:', df5c.shape)
            count_c, count_all = 0, 0
            for i in range(len(df_cells.columns)):
                c1 = df_cells.columns[i]
                count_c += 1
                c_list, p_list = [], []
                for j in range(len(df_drugs.columns)):
                    c2 = df_drugs.columns[j]
                    count_all += 1
                    s2 = df_drugs[c2] if argv.top in (0, None) else s_top(df_drugs[c2], argv.top)
                    ftols = pd.ols(y = df_cells[c1], x = s2, intercept=True)
                    c_list.append(df_cells[c1].corr(s2))
                    p_list.append(ftols.f_stat['p-value'])
                df5c[c1] = pd.Series(c_list).values
                df5p[c1] = np.array(p_list)
            print('Cells count: {}, Drugs count: {}'.format(count_c, count_all/count_c))
            # Output data frames to files
            df5c.to_csv('./corr_{}_{}.csv'.format(cells_base, drugs_base))
            df5p.to_csv('./p_value_{}_{}.csv'.format(cells_base, drugs_base))
            print('result shape:', df5c.shape)
        print('Done.')

    except Usage as err:
        print(sys.stderr, err.msg)
        print(sys.stderr, "Terminated, for help use -h or --help")
        return 2

if __name__ == "__main__":
    sys.exit(main())