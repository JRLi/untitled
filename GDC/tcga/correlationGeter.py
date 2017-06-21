#!/usr/bin/env python
import pandas as pd
import numpy as np
import argparse
import os, sys, datetime


def args_parse():
    parser = argparse.ArgumentParser(description='need correlation and pvalue dirs')
    parser.add_argument('-m', '--mode', type=str, choices=['p', 'c'], default='p', help="p or c, default is p")
    parser.add_argument('-t', '--threshold', type=float, default=0.05, help="threshold of cor or p, default is 0.1")
    args = parser.parse_args()
    return args


def openDF(in_path, direct='f'):
    fpath, fname = os.path.split(in_path)
    fbase, fext = os.path.splitext(fname)
    df = pd.read_csv(in_path, index_col=0) if fext == '.csv' else pd.read_table(in_path, index_col=0)
    if direct == 't':
        df = df.transpose()
    return df, fbase


def locate(df_input, mode, threshold):
    if mode == 'p':
        ii = np.where(df_input.values <= threshold)
    else:
        ii = np.where((df_input.values >= threshold) | (df_input.values <= -threshold))
    return ii


def lincs_dict(path):
    lineCount = 0
    id2anno = {}
    with open(path) as inputFile:
        for line in inputFile:
            lineCount += 1
            if lineCount == 1:
                continue
            lf = line.rstrip().split('\t')
            drug_annotation = '\t'.join([lf[1], lf[7], lf[6], lf[4], lf[9]]) if not lf[8].startswith('ctl_') \
                else '\t'.join([lf[1], '-', '-', '-', '-'])
            id2anno[lf[0]] = drug_annotation
    return id2anno


def main(argv=None):
    if argv is None:
        argv = args_parse()
        print(argv)
        id2drug = lincs_dict('GSE70138_20151231_annotation.txt')
        list_p = os.listdir('pvalue')
        time_1 = datetime.datetime.now()
        with open('result_{}_{}'.format(argv.mode, argv.threshold), 'a') as outFile:
            outFile.write('file\tlincsID\tcellLine\tdrug\tdrugID\tum\th\tsnvGene\tgeneID\tchr\tlocus\tcorrelation\tp_value\n')
            for fileName_p in list_p:
                fileName_c = fileName_p.replace('P_value', 'Corr', 1)
                path_p = os.path.join('pvalue', fileName_p)
                path_c = os.path.join('correlation', fileName_c)
                df_p, p_base = openDF(path_p)
                df_c, c_base = openDF(path_c)
                ii = locate(df_p, argv.mode, argv.threshold) if argv.mode == 'p' else locate(df_c, argv.mode, argv.threshold)
                title = p_base.replace('P_value_Corr_', '')
                if len(ii[1]) > 0:
                    for r, c in zip(ii[0], ii[1]):
                        snv = df_p.columns[c].replace(':', '\t')
                        drug = id2drug.get(df_p.index[r])
                        outFile.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                            title, df_p.index[r], drug, snv, df_c.iloc[r, c], df_p.iloc[r, c]))
                else:
                    print('{} has no value match the threshold: {}: {}'.format(title, argv.mode, argv.threshold))
        time_2 = datetime.datetime.now()
        print('\t[info] All used time:', str(time_2 - time_1))

if __name__ == "__main__":
    sys.exit(main())
