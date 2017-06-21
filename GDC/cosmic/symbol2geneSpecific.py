#!/usr/bin/env python
import argparse
import sys
from collections import defaultdict

HGNC_anno = 'D:/Project/drs/hgnc_id_170512'
pre_file = 'D:/Project/drs/gdsc/cosmic/cosv80_com_exp_df.tsv'
aft_file = 'D:/Project/drs/gdsc/cosmic/cosv80_com_EID_exp_df.tsv'
check_dict = 'D:/Project/drs/hgncSym2ID'
def main(argv=None):
    if argv is None:
        sym2entid_dict = {}
        with open(HGNC_anno) as hgnc_in:
            for line in hgnc_in:
                if line.startswith('hgnc_id'):
                    continue
                lf = line.replace('\r\n', '').replace('\n', '').replace('\"', '').split('\t')
                if lf[1] != '' and lf[4] != '':
                    sym2entid_dict[lf[1]] = lf[4]
                    if lf[3] != '':
                        pf = lf[3].split('|')
                        for prev in pf:
                            if prev not in sym2entid_dict:
                                sym2entid_dict[prev] = lf[4]
                    if lf[2] != '':
                        af = lf[2].split('|')
                        for alias in af:
                            if alias not in sym2entid_dict:
                                sym2entid_dict[alias] = lf[4]

        print(len(sym2entid_dict))
        with open(pre_file) as in_file, open(aft_file, 'w') as out_file:
            no_id_count = 0
            for line in in_file:
                if line.startswith('\t'):
                    out_file.write(line)
                    continue
                lf = line.replace('\r\n', '').replace('\n', '').replace('\"', '').split('\t', maxsplit=1)
                eid = sym2entid_dict.get(lf[0])
                if eid is not None:
                    out_file.write('\t'.join([eid, lf[1]]) + '\n')
                else:
                    print(lf[0])
                    no_id_count += 1
            print(no_id_count)

        with open(check_dict, 'w') as out_file:
            for k,v in sym2entid_dict.items():
                out_file.write('\t'.join([k, v]) + '\n')

if __name__ == "__main__":
    sys.exit(main())