#!/usr/bin/env python3.6
"""
usage: ./ppa_getMod.py default_matrix json_file
use shell script to perform all jason files with one default table
"""
import sys
import pandas as pd
import numpy as np
import os
import json


def prepare_output_dir(output_dir):
    if os.path.exists(output_dir):
        pass
    else:
        os.mkdir(output_dir)


def open_df(path_in, direct='n'):
    df_path, df_name = os.path.split(path_in)
    df_base, df_ext = os.path.splitext(df_name)
    df = pd.read_csv(path_in, index_col=0) if df_ext == '.csv' else pd.read_table(path_in, index_col=0)
    return df if direct == 'n' else df.T, df_base


def parse_json_to_matrix(mod_file, df1, df2, f_base):
    with open(mod_file, 'r') as json_file, open('summary_of_ppa', 'a') as summary:
        sub_matrix_list = json.load(json_file)
        m_c = 0
        for m_dict in sub_matrix_list:
            m_c += 1
            rows1_array = np.array(m_dict['rows1']) - 1
            rows2_array = np.array(m_dict['rows2']) - 1
            cols_array = np.array(m_dict['columns']) - 1
            df3 = df1.iloc[rows1_array, cols_array]
            df4 = df2.iloc[rows2_array, cols_array]
            len_r1 = 1 if isinstance(m_dict['rows1'], int) else len(m_dict['rows1'])
            len_r2 = 1 if isinstance(m_dict['rows2'], int) else len(m_dict['rows2'])
            len_c = 1 if isinstance(m_dict['columns'], int) else len(m_dict['columns'])
            df3.to_csv('mod_ppa/{}/mod{}_pc_{}x{}.csv'.format(f_base, m_c, len_r1, len_c))
            df4.to_csv('mod_ppa/{}/mod{}_lnc_{}x{}.csv'.format(f_base, m_c, len_r2, len_c))
            summary.write('{}\tmod{}\t{}\t{}\t{}\n'.format(f_base, m_c, len_r1, len_r2, len_c))
        return m_c


def main(argv=None):
    if argv is None:
        try:
            os.remove('summary_of_ppa')
        except OSError:
            pass
        prepare_output_dir('mod_ppa')
        f_list = os.listdir(sys.argv[3])
        f_list.sort()
        df1, df1_b = open_df(sys.argv[1])
        df2, df2_b = open_df(sys.argv[2])
        print("File1: {}\nFile2: {}".format(df1_b, df2_b))
        mod_c, mod_all = 0, 0
        for mod_file in f_list:
            mod_c += 1
            m_path = os.path.join(sys.argv[3], mod_file)
            js_base, js_ext = os.path.splitext(mod_file)
            prepare_output_dir('mod_ppa/{}'.format(js_base))
            mod_ic = parse_json_to_matrix(m_path, df1, df2, js_base)
            mod_all += mod_ic
            print('{}\t{}'.format(mod_ic, js_base))
        print('mod_file_processed:\t{}\nall_mod_count:\t{}'.format(mod_c, mod_all))


if __name__ == "__main__":
    sys.exit(main())
