#!/usr/bin/env python
import os
import sys
import argparse
import pandas as pd
import numpy as np
output_dir = './Split_zScore/'

use_message = '''
    splitByProject.py -p project_list -a column2project_annotation -i input_data_frame
'''

def args_parse():
    parser = argparse.ArgumentParser(description=use_message)
    parser.add_argument('-p', '--project', help='Project list file, one line with one name')
    parser.add_argument('-a', '--annotate', help='name to project file, txt file, NO title and index')
    parser.add_argument('-i', '--inputDF', help='data frame before split')
    args = parser.parse_args()
    return args


def prepare_output_dir():
    print("prepare output dir")
    if os.path.exists(output_dir):
        pass
    else:
        os.mkdir(output_dir)


def openDF(in_path, direct = 'f'):
    fpath, fname = os.path.split(in_path)
    fbase, fext = os.path.splitext(fname)
    df = pd.read_csv(in_path, index_col=0) if fext == '.csv' else pd.read_table(in_path, index_col=0)
    if direct == 't':
        df = df.transpose()
    return df, fbase


def main(argv=None):
    if argv is None:
        argv = args_parse()
    print(argv)
    prepare_output_dir()
    dfi, dfi_base = openDF(argv.inputDF)
    dfa = pd.read_table(argv.annotate)
    len_count = 0
    with open(argv.project) as project_file:
        project_list = [x.replace('\n', '') for x in project_file]
        for pid in project_list:
            ii = np.where(dfa.values == pid)[1]
            len_count += len(ii)
            df2 = dfi.iloc[:, ii]
            df2.to_csv('{}{}_{}.csv'.format(output_dir, dfi_base, pid))
            print('{}:\t{}\t'.format(pid, len(ii)))
        print('all count:', len_count)

if __name__ == "__main__":
    sys.exit(main())