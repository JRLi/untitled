#!/usr/bin/env python
import os
import sys
import scipy.stats as st
import numpy as np
import pandas as pd
import argparse

use_message = '''
    To transfer the expression profile matrix to up and down matrix for base algorithm.
    Need Python3 and numpy, pandas, scipy.
    If "NO" need median_normalizing, assign '-m'.
'''
out_suffix = '_up_down.txt'


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


def args_parse():
    parser = argparse.ArgumentParser(description=use_message)
    parser.add_argument('-t', '--threshold', type=int, default=10, help='log p-value threshold, default is 10')
    parser.add_argument("-i", "--indexMean", action="store_true", help="mean value according to index")
    parser.add_argument('-d', '--direct', type=str, choices=['n', 'r'],
                        default='n', help='value direction, n: normal, r: reverse(i.e. cmap), default is n')
    parser.add_argument("-m", "--median", action="store_false", help="median_normalizing mode, is set, no use it")
    parser.add_argument("-z", "--zscore", action="store_false", help="z score transfer, is set, no use it")
    parser.add_argument('profile', nargs='+', help="drug expression profile, csv or txt file list separated by space")
    args = parser.parse_args()
    return args


def openDF(in_path):
    fpath, fname = os.path.split(in_path)
    fbase, fext = os.path.splitext(fname)
    df = pd.read_csv(in_path, index_col=0) if fext == '.csv' else pd.read_table(in_path, index_col=0)
    return df, fbase


def median_normalizing(df_input):
    df = df_input.copy()
    print('median normalization')
    return df.sub(df.median(axis=1), axis=0)


def z_transfer(df_input, mode):
    df = df_input.copy()
    print('z score transfer.')
    # mean and std track column from row index 0 to end, so axis use default (index 0).
    avg = df.mean()
    std = df.std()
    for k in range(len(df.columns)):
        df.iloc[:, k] = (df.iloc[:, k] - avg[k]) / std[k]
    if mode == 'n':
        return df
    else:
        return -df


def z_to_p_log_trim_split(df_input, z_threshold, mode):
    df= df_input.copy()
    if mode == 'up':
        for k in range(len(df.columns)):
            if (k + 1) % 500 == 0:
                print(k + 1)
            df.iloc[:, k][df.iloc[:, k] < 0] = 0
            df.iloc[:, k] = -np.log10(st.norm.cdf(-df.iloc[:, k]) * 2)
            df.iloc[:, k][df.iloc[:, k] > z_threshold] = z_threshold
        df = df.add_suffix('_up')
    else:
        for k in range(len(df.columns)):
            if (k + 1) % 500 == 0:
                print(k + 1)
            df.iloc[:, k][df.iloc[:, k] > 0] = 0
            df.iloc[:, k] = -np.log10(st.norm.cdf(df.iloc[:, k]) * 2)
            df.iloc[:, k][df.iloc[:, k] > z_threshold] = z_threshold
        df = df.add_suffix('_dn')
    return df


def rescale(df_input):
    df = df_input.copy()
    minv = df.values.min()
    maxv = df.values.max()
    df = (df - minv) / (maxv - minv)
    df.columns = df.columns.str.replace('\s+', '')
    return df


def df_mean_index(df_input):
    gp = df_input.groupby(df_input.index)
    print('mean group by index.')
    return gp.mean()


def main(argv=None):
    try:
        if argv is None:
            argv = args_parse()
            print(argv)
            file_list = argv.profile
            for profile in file_list:
                df1, fileBase = openDF(profile)
                print(df1.shape)
                print(df1.columns[0])
                print(df1.index[0])
                df1 = df_mean_index(df1) if argv.indexMean else df1
                print(argv.indexMean, df1.shape)
                df1 = median_normalizing(df1) if argv.median else df1
                # df1.to_csv('test.csv')
                print(argv.median, df1.shape)
                mt = 'MN' if argv.median else 'noNM'
                df1 = z_transfer(df1, argv.direct) if argv.zscore else df1
                dir = '' if argv.direct == 'n' else '_rev'
                dfup = z_to_p_log_trim_split(df1, argv.threshold, 'up')
                dfdn = z_to_p_log_trim_split(df1, argv.threshold, 'dn')
                dfupdn = pd.concat([dfup, dfdn], axis=1)
                dfupdn = rescale(dfupdn)
                dfupdn.to_csv('./{}_{}_t{}{}{}'.format(fileBase, mt, str(argv.threshold), dir, out_suffix), sep='\t')


    except Usage as err:
        print(sys.stderr, err.msg)
        print(sys.stderr, "Terminated, for help use -h or --help")
        return 2

if __name__ == '__main__':
    sys.exit(main())