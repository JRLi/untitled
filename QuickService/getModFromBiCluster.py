#!/usr/bin/env python3.5
"""
usage: ./getModFromBiCluster.py default_matrix json_file
use shell script to perform all jason files with one default table
"""
import sys
import pandas as pd
import numpy as np
import os
import json


def prepare_output_dir(output_dir='.'):
    print("prepare output dir")
    if os.path.exists(output_dir):
        pass
    else:
        os.mkdir(output_dir)


def parse_json_to_subMatrix(mod_file, df1, fbase):
    with open(mod_file, 'r') as json_file, open('summary', 'a') as summary:
        subMatrix_list = json.load(json_file)
        mod_index = 1
        summary.write(fbase + '\n')
        for mod_dict in subMatrix_list:
            rows_array = np.array(mod_dict['rows']) - 1
            cols_array = np.array(mod_dict['columns']) - 1
            df2 = df1.iloc[rows_array, cols_array]
            df2.to_csv('{}/{}_{}_{}x{}.csv'.format(fbase,fbase,str(mod_index),str(len(rows_array)),str(len(cols_array))))
            summary.write('mod{}: rows:\t{}\tcolumns:\t{}\n'.format(str(mod_index),str(len(rows_array)),str(len(cols_array))))
            mod_index += 1


def main(argv=None):
    if argv is None:
        dfpath, dfname = os.path.split(sys.argv[1])
        dfbase, dfext = os.path.splitext(dfname)
        jspath, jsname = os.path.split(sys.argv[2])
        jsbase, jsext = os.path.splitext(jsname)
        prepare_output_dir(jsbase)
        print(dfext)
        df1 = pd.read_table(sys.argv[1], sep=',', index_col=0) if dfext == '.csv' else pd.read_table(sys.argv[1], index_col=0)
        print(df1.shape)
        parse_json_to_subMatrix(sys.argv[2], df1, jsbase)

if __name__ == "__main__":
    sys.exit(main())
