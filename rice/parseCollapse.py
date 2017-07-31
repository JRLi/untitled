#!/usr/bin/env python3.6
import argparse
import os
import sys
collapse = 'out_collapse'
miRNA = 'forMiRNA'
degradome = 'forDeg'


def args_parse():
    parser = argparse.ArgumentParser(description='For parse fa to miRNA or degradome')
    parser.add_argument('-m', '--miRNA', nargs=2, type=int, default=[18, 24], help='miRNA length range threshold')
    parser.add_argument('-d', '--degrade', type=int, default=5, help='degradome minimum length threshold')
    args = parser.parse_args()
    return args


def seq_parser(fa_file, mirna_t_list, deg_t):
    fbase, fext = os.path.splitext(fa_file)
    with open(os.path.join(collapse, fa_file)) as col_fa, open(os.path.join(miRNA, fbase + '_m' + fext), 'w') as \
            mir, open(os.path.join(degradome, fbase + '_d' + fext), 'w') as deg:
        seq_count, mir_count, deg_count, n_count = 0, 0, 0, 0
        for line in col_fa:
            if line.startswith('>'):
                seq_count += 1
                id = line.rstrip()
                seq = next(col_fa).rstrip()
                if ('N' in seq) or ('n' in seq):
                    n_count += 1
                    continue
                if (len(seq) >= mirna_t_list[0]) and (len(seq) <= mirna_t_list[1]):
                    mir.write(id + '\n' + seq + '\n')
                    mir_count += 1
                if len(seq) >= deg_t:
                    deg_count += 1
                    deg.write(id + '\n' + seq + '\n')
        print('File: {}\tSequence: {}\tmiRNA_fa: {}\tdegradome_fa: {}\tWith_N: {}'.format(
            fa_file, seq_count, mir_count, deg_count, n_count))


def main(argv=None):
    if argv is None:
        argv = args_parse()
        print(argv)
        fileList = os.listdir(collapse)
        file_count = 0
        for fa in fileList:
            file_count += 1
            seq_parser(fa, argv.miRNA, argv.degrade)
        print('Parse file: {}'.format(file_count))

if __name__ == '__main__':
    sys.exit(main())