#!/usr/bin/env python
import argparse
import sys
import os
import datetime
import pandas as pd
import secrets
import numpy as np
use_message = '''
    Need Python3.6 and numpy, pandas; or Anaconda.
    Usage: baseAnalysis.py [-t, -c] data reg.
    reg is a regulation profile performed by profile2updn.py.
    Example: nohup python -u baseAnalysis.py -q a -p 1000 data.csv reg_updn.txt 
'''


def args_parse():
    parser = argparse.ArgumentParser(description=use_message)
    parser.add_argument('-p', '--perm', type=int, default=100, help='Times of permutation, default is 100')
    parser.add_argument('-q', '--qType', type=str, choices=['a', 'm'],
                        default='m', help='quantile normalization type, a: mean, m: median, default is m')
    parser.add_argument("-n", "--no_median", action="store_true", help="no median_normalizing mode, is set, no use MN")
    parser.add_argument('-d', '--direct', choices=['n', 't'],
                        default=['n', 'n'], nargs=2, help="n is normal, t is transpose; default is n n")
    parser.add_argument('pairs', nargs=2, help="cell line (1) and drug (2) expression profile")
    args = parser.parse_args()
    return args


def prepare_output_dir(output_dir):
    if os.path.exists(output_dir):
        pass
    else:
        os.mkdir(output_dir)


def open_df(in_path, direct='f'):
    fpath, fname = os.path.split(in_path)
    fbase, fext = os.path.splitext(fname)
    df = pd.read_csv(in_path, index_col=0) if fext == '.csv' else pd.read_table(in_path, index_col=0)
    if direct == 't':
        df = df.transpose()
    return df, fbase


def quantile_normalize(df_input, q_type):
    df = df_input.copy()
    # compute rank
    dic = {}
    for col in df:
        dic.update({col: sorted(df[col])})
    sorted_df = pd.DataFrame(dic)
    print('To quantile normalization.\nqType:', 'mean' if q_type == 'a' else 'median')
    rank = sorted_df.mean(axis=1).tolist() if q_type == 'a' else sorted_df.median(axis=1).tolist()
    # sort
    for col in df:
        t = np.searchsorted(np.sort(df[col]), df[col])
        df[col] = [rank[i] for i in t]
    return df


def median_normalizing(df_input):
    print('To median normalization')
    df = df_input.copy()
    return df.sub(df.median(axis=1), axis=0)


def calculate_diff(cur_exp, cur_reg, row_num):
    fg1 = cur_reg.mul(cur_exp, axis=0).abs()
    bg1 = pd.DataFrame((1 - cur_reg)).mul(cur_exp, axis=0).abs()
    for i in range(1, len(fg1.index)):
        fg1.iloc[i, :] = fg1.iloc[i, :] + fg1.iloc[i - 1, :]
        bg1.iloc[i, :] = bg1.iloc[i, :] + bg1.iloc[i - 1, :]
    for i in range(len(fg1.columns)):
        fg1.iloc[:, i] = fg1.iloc[:, i] / fg1.iloc[row_num - 1, i]
        bg1.iloc[:, i] = bg1.iloc[:, i] / bg1.iloc[row_num - 1, i]
    diff = fg1 - bg1
    return diff


def calculate_es(data, reg):
    es = pd.DataFrame(0, data.columns, reg.columns)
    print('To calculate ES')
    for k in range(len(data.columns)):
        cur_exp = data.iloc[:, k].sort_values(ascending=False)
        cur_reg = reg.loc[cur_exp.index]
        xx = calculate_diff(cur_exp, cur_reg, len(data.index))
        pos_es = xx.max()
        neg_es = xx.min()
        pos_es[pos_es < 0] = 0
        neg_es[neg_es > 0] = 0
        # c= pos_es.where(pos_es > neg_es.abs()).fillna(0.0) + neg_es.where(neg_es.abs() >= pos_es).fillna(0.0)
        # for i in range(len(pos_es)):
        #    es.iloc[k,i] = pos_es[i] if pos_es[i] > neg_es.abs()[i] else neg_es[i]
        es.iloc[k, :] = pos_es[pos_es > neg_es.abs()].append(neg_es[neg_es.abs() >= pos_es]).reindex(index=pos_es.index)
    return es


def permutation(data, reg, perm):
    print('To permutation:', perm)
    cur_reg = reg
    # pos_es = pd.DataFrame(0, np.arange(len(reg.columns)), np.arange(perm))
    pos_es = pd.DataFrame(0, reg.columns, np.arange(perm))
    neg_es = pd.DataFrame(0, reg.columns, np.arange(perm))
    # cur_reg = cur_reg.reindex(index=['Lypla1', 'Atp6v1h', 'Fam150a', 'St18', 'Pcmtd1', 'Adhfe1', '3110035E14Rik'])
    for k in range(perm):
        se = secrets.choice(range(len(data.columns)))  # Random number from data columns number
        cur_exp = data.iloc[:, se].sample(frac=1)
        cur_exp = np.array(cur_exp)
        xx = calculate_diff(cur_exp, cur_reg, len(data.index))
        pos_es.iloc[:, k] = xx.max()
        neg_es.iloc[:, k] = xx.min()
    return pos_es, neg_es


def normalize_es(es, pos_es, neg_es):
    print('To normalize')
    pavg = pos_es.mean(axis=1)
    navg = neg_es.mean(axis=1).abs()
    # pos_npes = pos_es.div(pavg, axis=0)
    # neg_npes = neg_es.div(navg, axis=0)
    for k in range(len(es.index)):
        pos = es.iloc[k, :].div(pavg, axis=0)
        neg = es.iloc[k, :].div(navg, axis=0)
        es.iloc[k, :][es.iloc[k, :] > 0] = pos
        es.iloc[k, :][es.iloc[k, :] <= 0] = neg
    return es


def final_es(df_input):
    df_up = df_input.iloc[:, 0:len(df_input.columns) // 2]
    df_dn = df_input.iloc[:, len(df_input.columns) // 2:len(df_input.columns)]
    df_all = df_up - df_dn.values
    df_all.columns = [c[:-6] for c in df_all.columns]
    return df_all


def main(argv=None):
    if argv is None:
        argv = args_parse()
        time_1 = datetime.datetime.now()
        print('\nStart time:', str(time_1))
        prepare_output_dir('es_out')

        data, data_base = open_df(argv.pairs[0], argv.direct[0])
        reg, reg_base = open_df(argv.pairs[1], argv.direct[1])
        data = quantile_normalize(data, argv.qType)
        time_2 = datetime.datetime.now()
        print('Finish quantile normalization:', str(time_2 - time_1))

        qn_suffix = '_qnMedian' if argv.qType == 'm' else '_qnMean'
        ixs = reg.index.intersection(data.index)    # index order is followed later df
        # ixs = [ix for ix in data.index if ix in reg.index]
        print('The numbers of intersection features between data and reg:', len(ixs))
        data = data.loc[ixs]
        reg = reg.loc[ixs]
        data = data if argv.no_median else median_normalizing(data)
        nm_suffix = '_noNM' if argv.no_median else '_NM'
        es = calculate_es(data, reg)
        time_3 = datetime.datetime.now()
        print('Finish ES calculation, prepare to permutation:', str(time_3 - time_2))

        pos_es, neg_es = permutation(data, reg, argv.perm)
        es = normalize_es(es, pos_es, neg_es)
        es = es.add_suffix('.ES')
        print('es shape:', es.shape)
        time_4 = datetime.datetime.now()
        print('After permutation:', str(time_4 - time_3))
        # es.to_csv('./{}_{}{}{}_perm{}_ES.csv'.format(data_base, reg_base, qn_suffix, nm_suffix, argv.perm),
        # na_rep='NA')
        fes = final_es(es)
        fes.to_csv('es_out/{}_{}{}{}_p{}_fES.csv'.format(data_base, reg_base, qn_suffix, nm_suffix, argv.perm),
                   na_rep='NA')
        time_5 = datetime.datetime.now()
        print('Finished time:', str(time_5))
        print('All used time:', str(time_5 - time_1))


if __name__ == "__main__":
    sys.exit(main())
