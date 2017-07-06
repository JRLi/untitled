#!/usr/bin/env python3
import pandas as pd
import numpy as np
import scipy.stats as st
import sys, os, datetime, argparse
use_message = '''
    Need Python3 and numpy, pandas, scipy; or Anaconda.
    Usage: python -u indTtestSNV.py [-d d, -m m] -s df_snv -c [df_correlation1 ....]
'''


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


def args_parse():
    parser = argparse.ArgumentParser(description=use_message)
    parser.add_argument('-d', '--direct', choices=['n', 't'], default='n',
                        help="corr file direction, n is normal(columns are samples), t is transpose; default is n")
    parser.add_argument('-m', '--min', type=int, default=3, help='minimum number of samples per condition')
    parser.add_argument('-s', '--snv', help="snv-patient binary profile, columns need to be samples")
    parser.add_argument('-c', '--corr', nargs='+',
                        help="drug-patient correlation profile list, if columns aren't samples, set -d to t")
    args = parser.parse_args()
    return args


def prepare_output_dir(output_dir):
    if os.path.exists(output_dir):
        pass
    else:
        os.mkdir(output_dir)


def openDF(in_path, direct='n'):
    fpath, fname = os.path.split(in_path)
    fbase, fext = os.path.splitext(fname)
    df = pd.read_csv(in_path, index_col=0) if fext == '.csv' else pd.read_table(in_path, index_col=0)
    print('{}: Transpose: {}'.format(fbase, direct))
    if direct == 't':
        df = df.transpose()
    return df, fbase


def scipy_ttest_ind(s1, s2):
    return st.ttest_ind(s1, s2)


def t_by_index_of_df(df_snv, df_drug_corr, min):
    dfp = pd.DataFrame(index=df_drug_corr.index)
    count_snv, count_all_snv, count_drug = 0, 0, 0
    ii1_len_list = []
    for i in range(len(df_snv.index)):
        count_all_snv += 1
        ii1 = np.where(df_snv.iloc[i, :] == 1)
        ii0 = np.where(df_snv.iloc[i, :] != 1)
        ii1_len_list.append(str(len(ii1[0])))
        p_value_list = []
        if (len(ii1[0]) < min) or (len(ii0[0]) < min):
            print('\t[Skip]:', df_snv.index[i], len(ii1[0]), len(ii0[0]))
            continue
        for j in range(len(df_drug_corr.index)):
            count_drug += 1
            ss1 = df_drug_corr.iloc[j, ii1[0]]
            ss0 = df_drug_corr.iloc[j, ii0[0]]
            t_result = scipy_ttest_ind(ss1, ss0)
            p_value_list.append(t_result[1])
        count_snv += 1
        dfp[df_snv.index[i]] = np.array(p_value_list)
    print('[All_snv]: {}, [processed_snv]: {}, [drugs]: {}'.format(count_all_snv, count_snv,count_drug / count_snv))
    return dfp, ii1_len_list


def main(argv=None):
    try:
        if argv is None:
            argv = args_parse()
        time_0 = datetime.datetime.now()
        print('[Start]:{}\nsnv_file: {}\nmin: {}\ncorr_list: {}'.format(str(time_0), argv.snv, argv.min, argv.corr))
        prepare_output_dir('./p_value_df')

        df_snv, snv_base = openDF(argv.snv)
        print('{}: {}'.format(snv_base, df_snv.shape))

        with open('./Summary_snv_{}_{}'.format(snv_base, argv.min), 'w') as status:
            status.write('corr_file\tintersection\t{}\n'.format('\t'.join(df_snv.index)))
            for corr_file in argv.corr:
                time_1 = datetime.datetime.now()
                df_drug_corr, dc_base = openDF(corr_file, argv.direct)
                print('{}: {}\nExpect target shape: ({}, {})'.
                      format(dc_base, df_drug_corr.shape, df_drug_corr.shape[0], df_snv.shape[0]))

                ixc = df_snv.columns.intersection(df_drug_corr.columns)
                print('intersection samples:', len(ixc))
                df_corr = df_drug_corr[ixc]
                df_snv = df_snv[ixc]

                df_p, snv1_list = t_by_index_of_df(df_snv, df_corr, argv.min)
                print('Actual result shape:', df_p.shape)
                df_p.to_csv('p_value_df/ttest_{}_{}_{}.csv'.format(snv_base, dc_base, argv.min))
                status.write('{}\t{}\t{}\n'.format(dc_base, len(ixc), '\t'.join(snv1_list)))
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