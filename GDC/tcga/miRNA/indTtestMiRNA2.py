#!/usr/bin/env python3
import pandas as pd
import numpy as np
import scipy.stats as st
import sys
import os
import datetime
import argparse
use_message = '''
    Need Python3 and numpy, pandas, scipy; or Anaconda.
    Usage: python -u indTtestFormiRNA.py [-d d, -m m] -s df_label -c [df_miRNA_exp]
'''


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


def args_parse():
    parser = argparse.ArgumentParser(description=use_message)
    parser.add_argument('-d', '--direct', choices=['n', 't'], default='n',
                        help="exp file direction, n is normal(columns are samples), t is transpose; default is n")
    parser.add_argument('-n', '--nan', help='set specific value to np.nan')
    parser.add_argument('-v', '--variance', action="store_true", help='If set, perform a standard independent 2 sample '
                                                                      'test that assumes equal population variances')
    parser.add_argument('-f', '--fillna', action="store_true", help='set -f if need mean value to fill nan')
    parser.add_argument('-m', '--min', type=int, default=3, help='minimum number of samples per condition')
    parser.add_argument('-l', '--label', help="label-patient binary profile")
    parser.add_argument('-e', '--exp', nargs='+',
                        help="mirna-patient exp profile list, if columns aren't samples, set -d to t")
    args = parser.parse_args()
    return args


def prepare_output_dir(output_dir):
    if os.path.exists(output_dir):
        pass
    else:
        os.mkdir(output_dir)


def open_df(in_path, direct='n'):
    fpath, fname = os.path.split(in_path)
    fbase, fext = os.path.splitext(fname)
    df = pd.read_csv(in_path, index_col=0) if fext == '.csv' else pd.read_table(in_path, index_col=0)
    print('{}: Transpose: {}'.format(fbase, direct))
    if direct == 't':
        df = df.transpose()
    return df, fbase


def scipy_ttest_ind(s1, s2, var):
    return st.ttest_ind(s1, s2, equal_var=var)


def rm_zero(df_input):
    df = df_input.loc[(df_input != 0).any(1), (df_input != 0).any(0)]
    return df


def imputation_mean(df_input, direct=1):
    df = df_input.copy()
    df = df.T.fillna(df.T.mean()).T if direct == 1 else df.fillna(df.mean())
    return df


def t_by_index_of_df(df_label, df_exp, min_s, variance):
    dfp = pd.DataFrame(index=df_exp.index)
    count_mir = 0
    ii1_len_list = []
    ii1 = np.where(df_label.iloc[0, :] == 1)
    ii0 = np.where(df_label.iloc[0, :] == 0)
    ii1_len_list.append(str(len(ii1[0])))
    p_value_list, t_list_m, n_list_m, t_list_std, n_list_std= [], [], [], [], []
    if (len(ii1[0]) < min_s) or (len(ii0[0]) < min_s):
        print('\t[Skip]:', df_label.index[0], len(ii1[0]), len(ii0[0]))
        sys.exit(1)
    for i in range(len(df_exp.index)):
        count_mir += 1
        ss1 = df_exp.iloc[i, ii1[0]]
        ss0 = df_exp.iloc[i, ii0[0]]
        t_result = scipy_ttest_ind(ss1, ss0, variance)
        p_value_list.append(t_result[1])
        t_list_m.append(ss1.mean())
        n_list_m.append(ss0.mean())
        t_list_std.append(ss1.std())
        n_list_std.append(ss0.std())
    dfp['p_value'] = np.array(p_value_list)
    dfp['T_mean'] = np.array(t_list_m)
    dfp['N_mean'] = np.array(n_list_m)
    dfp['T_std'] = np.array(t_list_std)
    dfp['N_std'] = np.array(n_list_std)
    print('mir: {}'.format(count_mir))
    return dfp, ii1_len_list


def main(argv=None):
    try:
        if argv is None:
            argv = args_parse()
            print(argv)
        time_0 = datetime.datetime.now()
        print('[Start]:{}\nlabel_file: {}\nmin: {}\nexp_list: {}'.format(str(time_0), argv.label, argv.min, argv.exp))
        prepare_output_dir('./p_value_df')

        df_label, label_base = open_df(argv.label, 't')    # get label data frame
        print('{}: {}'.format(label_base, df_label.shape))

        with open('./Summary_{}_{}'.format(label_base, argv.min), 'w') as status:
            status.write('exp_file\tintersection\t{}\n'.format('\t'.join(df_label.index)))
            for exp_file in argv.exp:
                time_1 = datetime.datetime.now()
                df_exp, dc_base = open_df(exp_file, argv.direct)     # get exp data frame

                if argv.nan:
                    df_exp = df_exp.replace(argv.nan, np.nan)
                    df_exp = df_exp.astype(float)

                if argv.fillna:
                    print('[Process fill nan, default]:', df_exp.isnull().sum(1).sum(), sep='\n')
                    df_exp = imputation_mean(df_exp, 1) if argv.direct == 'n' else imputation_mean(df_exp)
                    print('[Process fill nan, after]:', df_exp.isnull().sum(1).sum(), sep='\n')
                print('{}: {}\nExpect target shape: ({}, {})'.
                      format(dc_base, df_exp.shape, df_exp.shape[0], df_label.shape[0]))

                ixc = df_label.columns.intersection(df_exp.columns)     # columns intersection
                print('intersection samples:', len(ixc))
                df_exp = df_exp[ixc]
                df_label = df_label[ixc]

                df_p, label_1_list = t_by_index_of_df(df_label, df_exp, argv.min, argv.variance)
                suffix = 'standard_ind_t' if argv.variance else "welch_t"
                print('Actual result shape:', df_p.shape)
                df_p.to_csv('p_value_df/ttest_{}_{}_{}.csv'.format(dc_base, argv.min, suffix))
                status.write('{}\t{}\t{}\n'.format(dc_base, len(ixc), '\t'.join(label_1_list)))

                time_3 = datetime.datetime.now()
                print('\t[Finished time]: {}\t[Used time]: {}'.format(str(time_3), str(time_3 - time_1)))

            time_4 = datetime.datetime.now()
            print('\t[All finish time]: {}\t[All use time]: {}\nDone.\n'.format(str(time_4), str(time_4 - time_0)))

    except Usage as err:
        print(sys.stderr, err.msg)
        print(sys.stderr, "Terminated, for help use -h or --help")
        return 2

if __name__ == '__main__':
    sys.exit(main())
