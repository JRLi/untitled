#!/usr/bin/env python3
import os
import pandas as pd


def open_df(in_path, direct='n'):
    fpath, fname = os.path.split(in_path)
    fbase, fext = os.path.splitext(fname)
    df = pd.read_csv(in_path, index_col=0) if fext == '.csv' else pd.read_table(in_path, index_col=0)
    print('{}: Transpose: {}'.format(fbase, direct))
    if direct == 't':
        df = df.transpose()
    return df, fbase


def get_3list(word):
    p_list, c_list, l_list = [], [], []
    for proj in word.split(' '):
        print(proj)
        with open(proj.join(['label/lab_TCGA-', '.csv'])) as in_f:
            next(in_f)
            for line in in_f:
                case, check = line.rstrip().split(',')
                label = 'Cancer' if check == '1' else 'Normal'
                p_list.append('TCGA-' + proj)
                c_list.append(case)
                l_list.append(label)
    return p_list, c_list, l_list

def main():
    df1, df_base = open_df('./mirna_rpm_caseID.csv')
    #word1 = 'BRCA KIRC THCA PRAD LIHC LUAD LUSC STAD HNSC KIRP UCEC KICH BLCA ESCA COAD'
    word1 = 'COAD'
    proj_list, case_list, label_list = get_3list(word1)
    print(len(proj_list), len(case_list), len(label_list))
    df1 = df1[case_list]
    df1.to_csv(df_base + '_15.csv')
    with open('first_two_row.csv', 'w') as out_f:
        out_f.write('Cancer,' + ','.join(proj_list) + '\n' + 'Type,' + ','.join(label_list) + '\n')

if __name__ == '__main__':
    main()