#!/usr/bin/env python3
from collections import defaultdict
import os, sys


def data_files_dict(in_path):
    d2f = defaultdict(list)
    f2d = defaultdict(set)
    with open(in_path) as in_file:
        for line in in_file:
            if line.startswith('D00'):
                lf = line.rstrip().split('\t')
                d2f[lf[0]].append(lf[1])
                f2d[lf[1]].add(lf[0])
    return d2f, f2d


def main():
    data2files, file2data = data_files_dict('ImmSamples.tsv')
    print('find common samples in ImmSamples.tsv, all samples:', len(file2data))
    for k, v in file2data.items():
        c2 =0
        if len(v) > 1:
            c2 += 1
            print(k, v, sep='\t')
    print('{} samples are shared in different data sets')
    print('check downloaded files')
    dir_count = 0
    for dirPath, dirNames, fileNames in os.walk('.'):
        if dirPath.startswith('./D00'):
            dir_count += 1
            file_list = [x.replace('.tsv', '') for x in fileNames]
            data = dirPath.replace('./', '')
            sample_list = data2files.get(data)
            print('{}: should have samples: {}\tactual files:{}'.format(data, len(sample_list), len(file_list)))
            if len(sample_list) != len(file_list):
                print('\t[Check]set of samples:{}\tset of files:{}'.format(len(set(sample_list)), len(set(file_list))))
            print('symmetric_difference:', set(sample_list).symmetric_difference(file_list))
    print('all data:{}'.format(dir_count))

if __name__ == '__main__':
    main()
