#!/usr/bin/env python
import sys
import os
import numpy as np
import pandas as pd
import argparse
import scipy.stats as st

use_message = '''
    Need Python3 and numpy, pandas, scipy; or Anaconda.
    Usage: corrV3.py [-t, -c] df1 df2
    Example: corrV3.py -t 250 -c s df1.csv df2.txt 
    OR: corrV3.py df1 df2, will use all genes and pearson correlation.
'''
mRNA_d = 'Split_mRNA'
miRNA_d = 'split'
mRNA_lab_d = 'label_mRNA'
miRNA_lab_d = 'label'


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


def args_parse():
    parser = argparse.ArgumentParser(description=use_message)
    parser.add_argument('-c', '--corr', type=str, choices=['p', 'k', 's'],
                        default='p', help='Correlation method, p: pearson, k: kendall, s: spearman; default is p')
    parser.add_argument('-d', '--direct', choices=['n', 't'],
                        default=['n', 'n'], nargs=2, help="n is normal, t is transpose; default is n n")
    parser.add_argument('-t', '--type', type=str, choices=['t', 'n', 'a'], default='t',
                        help='type of samples, t: tumor, n: normal, a: all, default is t')
    parser.add_argument('pairs', nargs=2, help="mRNA (1) and miRNA (2)")
    args = parser.parse_args()
    return args


def open_df(in_path, direct):
    fpath, fname = os.path.split(in_path)
    fbase, fext = os.path.splitext(fname)
    df = pd.read_csv(in_path, index_col=0) if fext == '.csv' else pd.read_table(in_path, index_col=0)
    print('{}: Transpose: {}'.format(fbase, direct))
    if direct == 't':
        df = df.transpose()
    return df, fbase


def scipy_corr(s1, s2, corr_mode):
    ixs = s1.index.intersection(s2.index)
    if corr_mode == 'pearson':
        return st.pearsonr(s1[ixs], s2[ixs])
    elif corr_mode == 'kendall':
        return st.kendalltau(s1[ixs], s2[ixs])
    elif corr_mode == 'spearman':
        return st.spearmanr(s1[ixs], s2[ixs])


def main(argv=None):
    try:
        if argv is None:
            argv = args_parse()
            mRNAs = os.listdir(mRNA_d)
            tcga_count = 0
            final_dict = {}
            col_list = ['correlation', 'p_value', 'mRNA_tumor', 'miRNA_tumor', 'intersection']
            s_type = {'t': 'tumor', 'n': 'normal', 'a': 'all'}.get(argv.type)
            for m_file in mRNAs:
                if m_file.startswith('TCGA'):
                    tcga_count += 1
                    tmp_dict = {}

                    dfm, dfm_base = open_df(os.path.join(mRNA_d, m_file), argv.direct[0])
                    dfmi, dfmi_base = open_df(os.path.join(miRNA_d, dfm_base + '_mirna_rpm_caseID.csv'), argv.direct[1])
                    dfmi.columns = ['-'.join(x.split('-')[0:4]) for x in dfmi.columns]
                    #dfmi.to_csv('test.csv')
                    dfml, lm_base = open_df(os.path.join(mRNA_lab_d,'lab_' + dfm_base + '.csv'), 't')
                    dfmil, lmi_base = open_df(os.path.join(miRNA_lab_d,'lab_' + dfm_base + '.csv'), 't')

                    if  s_type == 'all':
                        s_m = dfm.ix[int(argv.pairs[0]), :]
                        s_mi = dfmi.ix[argv.pairs[1], :]
                    elif s_type == 'tumor':
                        ii1_m = np.where(dfml.iloc[0, :] == 1)
                        ii1_mi = np.where(dfmil.iloc[0, :] == 1)
                        s_m = dfm.ix[int(argv.pairs[0]), ii1_m[0]]
                        s_mi = dfmi.ix[argv.pairs[1], ii1_mi[0]]
                    else:
                        ii1_m = np.where(dfml.iloc[0, :] == 0)
                        ii1_mi = np.where(dfmil.iloc[0, :] == 0)
                        s_m = dfm.ix[int(argv.pairs[0]), ii1_m[0]]
                        s_mi = dfmi.ix[argv.pairs[1], ii1_mi[0]]
                    mt_num = len(s_m)
                    mit_num = len(s_mi)
                    print('[{}]\tmRNA_tumor:{}\tmiRNA_tumor:{}'.format(dfm_base, mt_num, mit_num))

                    s_mi = s_mi[~s_mi.index.duplicated(keep='first')]   # super important, remove duplicated index
                    ixs = s_m.index.intersection(s_mi.index)
                    print('The numbers of intersection samples:', len(ixs))
                    s_m = s_m[ixs]
                    s_mi = s_mi[ixs]

                    if dfm_base == 'TCGA-BRCA':
                        df1 = pd.DataFrame()
                        df1['g1028'] = s_m
                        df1[argv.pairs[1]] = s_mi
                        df1.to_csv('BRCA_{}_{}_{}.csv'.format(argv.pairs[1], 'g1028', s_type))

                    method = {'p': 'pearson', 'k': 'kendall', 's': 'spearman'}.get(argv.corr)
                    cor_e, p_value = scipy_corr(s_m, s_mi, method)
                    tmp_dict[col_list[0]] = cor_e
                    tmp_dict[col_list[1]] = p_value
                    tmp_dict[col_list[2]] = mt_num
                    tmp_dict[col_list[3]] = mit_num
                    tmp_dict[col_list[4]] = len(ixs)
                    final_dict[dfm_base] = tmp_dict
            df = pd.DataFrame(final_dict)
            df = df.reindex(index=col_list)
            df.T.to_csv('{}_{}_{}_correlation.csv'.format(argv.pairs[1], argv.pairs[0], s_type))
            print('Done. process {} files.\n'.format(tcga_count))

    except Usage as err:
        print(sys.stderr, err.msg)
        print(sys.stderr, "Terminated, for help use -h or --help")
        return 2

if __name__ == "__main__":
    sys.exit(main())