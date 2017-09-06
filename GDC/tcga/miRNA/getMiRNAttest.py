#!/usr/bin/env python3
import os
import sys
import argparse


def args_parse():
    parser = argparse.ArgumentParser(description='need pvalue dir')
    parser.add_argument('-l', '--level', type=str, choices=['r', 'c'], default='r', help='rpm(r), or count(c)')
    parser.add_argument('-m', '--miRNAs', type=str, nargs='+', help='miRNA name, separate by space')
    parser.add_argument('-s', '--status', help='summary file of miRNA label, e.g. check_gdc')
    args = parser.parse_args()
    return args


def p2cn_dict(path):
    pro2cn_dict = {}
    with open(path) as in_f:
        next(in_f)
        for line in in_f:
            lf = line.rstrip().split('\t', maxsplit=2)
            pro2cn_dict[lf[0]] = lf[2].replace('\t', ',')
        return  pro2cn_dict


def main(argv=None):
    if argv is None:
        argv = args_parse()
        pro2cn_dict = p2cn_dict(argv.status)
        p_folder = 'p_value_df' if argv.level == 'r' else 'p_value_df_rc'
        list_p = os.listdir(p_folder)

        with open('miRNA_extraction.csv', 'w') as out_f:
            out_f.write('miRNA,Cancer_type,P_value,mean_rpm_tumor,mean_rpm_normal,std_rpm_tumor,std_rpm_normal,#tumor_case,#normal_case\n')
            file_count = 0
            for fileName_p in list_p:
                if not fileName_p.endswith('_welch_t.csv'):
                    continue
                file_count += 1
                check_count = 0
                p_id = fileName_p.split('_')[1]
                tn = pro2cn_dict.get(p_id)
                with open(os.path.join(p_folder, fileName_p)) as in_f:
                    next(in_f)
                    for line in in_f:
                        lf = line.rstrip().split(',')
                        if lf[0] in argv.miRNAs:
                            check_count += 1
                            out_f.write('{},{},{},{},{},{},{},{}\n'.format(lf[0], p_id, lf[1], lf[2], lf[3], lf[4], lf[5], tn))
                    print('[{}] check_miRNAs:{}\tchecked: {}'.format(p_id, len(argv.miRNAs),check_count))
            print('Done: {} files'.format(file_count))

if __name__ == '__main__':
    sys.exit(main())
