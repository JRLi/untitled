#!/usr/bin/env python3.6
import os
import pandas as pd
from collections import defaultdict


def get_dict_list(path, file_n):
    m2r_dict = defaultdict(list)
    read_list = []
    with open(os.path.join(path, file_n)) as in_f:
        next(in_f)
        for line in in_f:
            lf = line.rstrip().split('\t')
            read_list.append(lf[1])
            m2r_dict[lf[3]].append(lf[1])
        print('\t{}: {} miRNA'.format(file_n, len(m2r_dict)))
        return read_list, dict(sorted(m2r_dict.items()))


def process_mir2exp_dict(ar_list, m2rl_dict):
    mir2exp_dict = {}
    for miRNA, r_list in m2rl_dict.items():
        exp = 0
        for read_n in r_list:
            count_r = int(read_n.split('_x')[1])
            exp += count_r / ar_list.count(read_n)
        mir2exp_dict[miRNA] = exp
    return mir2exp_dict


def get_df(path, f_list, key_wd):
    file2m2e_dict = {}     # key is file name, value is miRNA to expression dictionary
    ct_file = 0
    for file_n in f_list:
        if key_wd in file_n:    # process files of check
            ct_file += 1
            fbase, fext = os.path.splitext(file_n)
            ar_list, m2rl_dict = get_dict_list(path, file_n)  # get dict and list for single file
            m2exp_dict = process_mir2exp_dict(ar_list, m2rl_dict)
            file2m2e_dict[fbase] = m2exp_dict
    print('[Info]{}: {} files'.format(key_wd, ct_file))
    return file2m2e_dict, ct_file


def main():
    check_list = ['CD4', 'CD8', 'CD14', 'CD15', 'CD19', 'CD56', 'CD235a', 'WB', 'serum', 'exosome']
    path_d = 'GSE100467'
    file_list = os.listdir(path_d)
    f_count = 0
    for check in check_list:    # 10 check
        print('[Start process]', check)
        c2m2e_dict, file_count = get_df(path_d, file_list, check)
        df1 = pd.DataFrame(c2m2e_dict)
        df1 = df1.fillna(0)
        df1.to_csv(path_d.lower() + '_' + check + '.csv')
        f_count += file_count
    print('[Done], all files:', f_count)

if __name__ == '__main__':
    main()