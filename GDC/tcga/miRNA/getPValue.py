#!/usr/bin/env python3.6
import os
import sys
import argparse


def args_parse():
    parser = argparse.ArgumentParser(description='need pvalue dir')
    parser.add_argument('-t', '--threshold', type=float, default=0.01, help="threshold of p-value, default is 0.1")
    parser.add_argument('-l', '--level', type=str, choices=['r', 'c'], default='r', help='rpm(r), or count(c)')
    parser.add_argument('-s', '--status', help='summary file of miRNA label')
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

        with open('p_value_{}_results.csv'.format(argv.threshold), 'w') as out_f:
            out_f.write('Cancer_type,miRNA,#tumor_case,#normal_case,mean_rpm_tumor,mean_rpm_normal,P_value\n')
            file_count, all_check = 0, 0
            for fileName_p in list_p:
                if not fileName_p.startswith('ttest'):
                    continue
                file_count += 1
                line_count, check_count = 0, 0
                p_id = fileName_p.split('_')[1]
                tn = pro2cn_dict.get(p_id)
                path_p = os.path.join(p_folder, fileName_p)
                with open(path_p) as in_f:
                    next(in_f)
                    for line in in_f:
                        line_count += 1
                        lf = line.rstrip().split(',')
                        if float(lf[1]) <= argv.threshold:
                            check_count += 1
                            out_f.write('{},{},{},{},{},{}\n'.format(p_id, lf[0], tn, lf[2], lf[3], lf[1]))
                    all_check += check_count
                    print('[{}] miRNA: {}\tcheck: {}'.format(p_id, line_count, check_count))
            print('Done: {} files, get {} miRNA'.format(file_count, all_check))

if __name__ == '__main__':
    sys.exit(main())