#!/usr/bin/env python
import pandas as pd
import numpy as np
import argparse
import os, sys, datetime


def args_parse():
    parser = argparse.ArgumentParser(description='need pvalue dir')
    parser.add_argument('-t', '--threshold', type=float, default=0.05, help="threshold of p-value, default is 0.1")
    parser.add_argument('-l', '--level', type=str, choices=['g', 'l'], default='g', help='snv level, default is g')
    parser.add_argument('-s', '--status', help='summary file output by t-test')
    parser.add_argument('-p', '--p_value', default='p_value_df', help='p_value directory')
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


def snv2number_dict(path):
    lineCount, intersect = 0, 0
    snv_list, num_list = [], []
    with open(path) as inputFile:
        for line in inputFile:
            lf = line.rstrip().split('\t')
            lineCount += 1
            if lineCount == 1:
                snv_list = lf[2:]
            elif lineCount == 2:
                intersect = lf[1]
                num_list = lf[2:]
            else:
                break
        snv2number = {k: v for k, v in zip(snv_list, num_list)}
    return snv2number, intersect


def rreplace(s, old, new, occurrence):
    li = s.rsplit(old, occurrence)
    return new.join(li)


def main(argv=None):
    if argv is None:
        argv = args_parse()
        print(argv)
        id2drug = lincs_dict('GSE70138_20151231_annotation.txt')
        list_p = os.listdir(argv.p_value)
        time_1 = datetime.datetime.now()
        snv2num, int_p = snv2number_dict(argv.status)

        with open('pvalue_exp_vs_lincs_vs_snv_{}_{}.txt'.format(argv.level, argv.threshold), 'w') as outFile:
            outFile.write('corr_file\tGDC\tlincsID\tcellLine\tdrug\tdrugID\tum\th\tsnvGene\tgeneID\t')
            outFile.write('chr\tlocus\tinter_p\tsnv_num\tp_value\n' if argv.level == 'l' else 'inter_p\tsnv_num\tp_value\n')
            for fileName_p in list_p:
                path_p = os.path.join(argv.p_value, fileName_p)
                df_p, p_base = openDF(path_p)
                title = p_base.replace('snv_indel.extracted.', '')
                title = rreplace(title, '.mut.status', '', 1)
                gdc = title[title.find('snp.') + 4:]
                gdc = gdc[:gdc.find('.')]

                ii = locate(df_p, 'p', argv.threshold)
                if len(ii[1]) > 0:
                    for r, c in zip(ii[0], ii[1]):
                        snv_num = snv2num.get(df_p.columns[c])
                        snv = df_p.columns[c].replace(':', '\t')
                        drug = id2drug.get(df_p.index[r])
                        outFile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                            title, gdc, df_p.index[r], drug, snv, int_p, snv_num, df_p.iloc[r, c]))
                else:
                    print('{} has no value match the threshold: {}: {}'.format(title, argv.mode, argv.threshold))
        time_2 = datetime.datetime.now()
        print('\t[info] All used time:', str(time_2 - time_1))

if __name__ == "__main__":
    sys.exit(main())
