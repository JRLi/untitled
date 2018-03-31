#!/usr/bin/env python3.6
import os
import sys
import argparse
import pandas as pd
use_message = ''


def args_parse():
    parser = argparse.ArgumentParser(description=use_message)
    parser.add_argument('-a', '--axis', choices=['c', 'i'], default='c', help='The axis to concatenate along, c is '
                                                                              'columns/1, i is index/0, default is c')
    parser.add_argument('-j', '--join', choices=['o', 'i'], default='o', help='concat join strategy, o is outer, '
                                                                              'i is inner, default is o')
    parser.add_argument('-o', '--output', type=str, default='concat_out', help='output prefix, default is concat_out')
    parser.add_argument('files', nargs='+', help="concat data frame")
    args = parser.parse_args()
    return args


def open_df(in_path, direct):
    fpath, fname = os.path.split(in_path)
    fbase, fext = os.path.splitext(fname)
    df = pd.read_csv(in_path, index_col=0) if fext == '.csv' else pd.read_table(in_path, index_col=0)
    if direct == 't':
        df = df.transpose()
    return df, fbase


def df_stat(df_input, direct='col_wise'):
    df = df_input.copy()
    df2 = pd.DataFrame()
    axis_dir = 0 if direct == 'col_wise' else 1
    df2['avg'] = df.mean(axis_dir)
    df2['std'] = df.std(axis_dir)
    df2['max'] = df.max(axis_dir)
    df2['min'] = df.min(axis_dir)
    return df2


def df_concat(files, axis, join):
    ax = 1 if axis == 'c' else 0
    jn = 'outer' if join == 'o' else 'inner'
    appended_data = []
    c_count = 0
    for in_file in files:
        df, df_b = open_df(in_file, 'n')
        print('{}:{}'.format(df.shape, df_b))
        c_count += df.shape[ax]
        appended_data.append(df)
    print('all concat: {}'.format(c_count))
    appended_data = pd.concat(appended_data, axis=ax, join=jn)
    return appended_data


def main(argv=None):
    if argv is None:
        argv = args_parse()
        df = df_concat(sorted(argv.files), argv.axis, argv.join)
        print('result shape: {}'.format(df.shape))
        df.to_csv(argv.output + '.csv')


if __name__ == '__main__':
    sys.exit(main())
