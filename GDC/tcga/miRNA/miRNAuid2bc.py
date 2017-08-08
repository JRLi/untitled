#!/usr/bin/env python3.6
import argparse
import os
import sys
import pandas as pd
from collections import defaultdict


def argv_parse():
    parser = argparse.ArgumentParser(description='uuid -a anno file1 file2 .....')
    parser.add_argument('-a', '--annotation', help='annotation file')
    parser.add_argument('mirDataFrame', nargs='+', help='miRNA data frame files')
    args = parser.parse_args()
    return args


def prepare_output_dir(dir):
    print("prepare output dir")
    if os.path.exists(dir):
        pass
    else:
        os.mkdir(dir)


def openDF(in_path, direct = 'f'):
    fpath, fname = os.path.split(in_path)
    fbase, fext = os.path.splitext(fname)
    df = pd.read_csv(in_path, index_col=0) if fext == '.csv' else pd.read_table(in_path, index_col=0)
    if direct == 't':
        df = df.transpose()
    return df, fbase


def uuid2case_dict(path):
    with open(path) as in_f:
        uuid2cases, pro2case_l, l_count = {}, defaultdict(list), 0
        next(in_f)
        for line in in_f:
            l_count += 1
            lf = line.rstrip().split(',')
            uuid2cases[lf[0]] = lf[2]
            pro2case_l[lf[3]].append(lf[2])
        print('l_count: {}\tlen(uuid2cases): {}\tlen(pro2case): {}'.format(l_count, len(uuid2cases), len(pro2case_l)))
        return uuid2cases, pro2case_l


def split_df(df_in, p2cl_dict, d_base):
    print('Start to split files')
    for project, case_list in p2cl_dict.items():
        df = df_in[case_list]
        print('[info]{}: {}'.format(project, df.shape))
        df.to_csv('split/{}_{}_caseID.csv'.format(project, d_base))
        df.T.to_csv('transpose_split/{}_{}_caseID.csv'.format(project, d_base))


def main(argv=None):
    if argv is None:
        argv = argv_parse()
        uid2case_dict, pro2case_l_dict = uuid2case_dict(argv.annotation)
        prepare_output_dir('split')
        prepare_output_dir('transpose_split')

        for f_path in argv.mirDataFrame:
            df1, df_base = openDF(f_path)
            c_list = []
            for i in range(len(df1.columns)):
                cID = uid2case_dict.get(df1.columns[i], 'unknown-case-id-xxx-NN')
                c_list.append(cID)
                if cID == 'unknown-case-id-xxx-0':
                    print(df1.columns[i], 'have no barcode')
            df1.columns = c_list
            df1.to_csv(df_base + '_caseID.csv')

            split_df(df1, pro2case_l_dict, df_base)

if __name__ == '__main__':
    sys.exit(main())