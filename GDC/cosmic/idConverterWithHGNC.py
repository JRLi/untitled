#!/usr/bin/env python
import argparse
import sys
from collections import defaultdict

HGNC_anno = 'D:/Project/drs/hgnc_id_170512'


def args_parse():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-p', '--perm', type=int, default=100, help='Times of permutation, default is 100')
    parser.add_argument('-q', '--qType', type=str, choices=['a', 'm'],
                        default='m', help='quantile normalization type, a: mean, m: median, default is m')
    parser.add_argument("-n", "--no_median", action="store_true", help="no median_normalizing mode, is set, no use MN")
    parser.add_argument('-d', '--direct', choices=['n', 't'],
                        default=['n', 'n'], nargs=2, help="n is normal, t is transpose; default is n n")
    parser.add_argument('pairs', nargs=2, help="cell line (1) and drug (2) expression profile")
    args = parser.parse_args()
    return args


def main(argv=None):
    if argv is None:
        sym2entid_dict = defaultdict(list)
        with open(HGNC_anno) as hgnc_in:
            for line in hgnc_in:
                if line.startswith('hgnc_id'):
                    continue
                lf = line.replace('\r\n', '').replace('\n', '').replace('\"', '').split('\t')
                if lf[1] != '' and lf[4] != '':
                    sym2entid_dict[lf[1]].append(lf[4])
                    if lf[2] != '':
                        af = lf[2].split('|')
                        for alias in af:
                            sym2entid_dict[alias].append(lf[4])
                    if lf[3] != '':
                        pf = lf[3].split('|')
                        for prev in pf:
                            sym2entid_dict[prev].append(lf[4])

        print(len(sym2entid_dict))
        ambiguous = 0
        for k, v in sym2entid_dict.items():
            if len(v) > 1:
                ambiguous += 1
                print(k, v)     # no result, indicate all symbol to gene_id is one-by-one
        print(ambiguous)
        with open('D:/Project/drs/gdsc/cosmic/cosv81_raw_exp_df.tsv') as in_file, \
                open('D:/Project/drs/gdsc/cosmic/cosv81_raw_EID_exp_df.tsv', 'w') as out_file:
            no_id_count = 0
            for line in in_file:
                if line.startswith('\t'):
                    out_file.write(line)
                    continue
                lf = line.replace('\r\n', '').replace('\n', '').replace('\"', '').split('\t', maxsplit=1)
                eid_list = sym2entid_dict.get(lf[0])
                if eid_list is not None:
                    for eid in eid_list:
                        out_file.write('\t'.join([eid, lf[1]]) + '\n')
                else:
                    print(lf[0])
                    no_id_count += 1
            print(no_id_count)

if __name__ == "__main__":
    sys.exit(main())