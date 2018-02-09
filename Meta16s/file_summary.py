#!/usr/bin/env python3.6
import os
import numpy as np


def sample_type_dict(f_path):
    s2dx1 = {}
    with open(f_path) as in_f:
        next(in_f)
        for line in in_f:
            lf = line.split('\t')
            s2dx1[lf[0].replace('sample', '')] = lf[4]
        return s2dx1


def file_summary(path):
    with open(path) as in_f:
        lc, rc, r_len, = 0, 0, 0
        rl_list = []
        for line in in_f:
            lc += 1
            if lc % 4 == 2:
                rc += 1
                rl = len(line.rstrip())
                r_len += rl
                rl_list.append(rl)
        rl_ar = np.array(rl_list)
        if rc != lc/4:
            print(path, 'have line count problem')
        return rc, np.mean(rl_ar), np.std(rl_ar)


def summary_dict(path):
    f_list = os.listdir(path)
    ss_dict = {}
    f_c = 0
    for f_n in f_list:
        if 'R1' in f_n:
            f_c += 1
            prefix = f_n[: f_n.find('.')]
            rc_r1, avg_len_r1, std_len_r1 = file_summary(os.path.join(path, f_n))
            rc_r2, avg_len_r2, std_len_r2 = file_summary(os.path.join(path, prefix + '.R2.clean.fastq'))
            rc_ra, avg_len_ra, std_len_ra = file_summary('../../assemblied_reads/' + prefix +'.assembled.fastq')
            if rc_r1 != rc_r2:
                print('{} have {} R1, but have {} R2'.format(prefix, rc_r1, rc_r2))
            s = '{}\t{}\t{}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}'.\
                format(rc_ra, rc_r1, rc_r2, avg_len_ra, avg_len_r1, avg_len_r2, std_len_ra, std_len_r1, std_len_r2)
            ss_dict[prefix] = s
    print('files: {}'.format(f_c))
    return ss_dict


def main():
    s2t = sample_type_dict('ra-s-metadata.tsv')
    s2s = summary_dict('.')
    with open('read_check', 'w') as out_f:
        out_f.write('ID\tDx1\tread_m\tread_r1\tread_r2\tlen_m\tlen_r1\tlen_r2\tlen_std_m\tlen_std_r1\tlen_std_r2\n')
        for k, v in s2s.items():
            out_f.write('{}\t{}\t{}\n'.format(k, s2t.get(k, 'NaN'), v))


if __name__ == '__main__':
    main()
