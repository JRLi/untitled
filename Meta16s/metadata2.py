#!/usr/bin/env python3.6
import os
import secrets
ra_dir = 'D:/Project/106ra'


def manifest_set(path_i):
    with open(path_i) as in_f:
        sample_set = set()
        next(in_f)
        for line in in_f:
            lf = line.rstrip().split(',')
            sample_set.add(lf[0])
        return sample_set


def modify_man(path_i, c_l):
    with open(path_i) as in_f, open(os.path.join(ra_dir, 'ra-assemb-manifest_182f'), 'w') as out_f:
        out_f.write(next(in_f))
        for line in in_f:
            lf = line.rstrip().split(',')
            if lf[0] in c_l:
                out_f.write(line)
            else:
                print(lf[0])


def main():
    s_set = manifest_set(os.path.join(ra_dir, 'ra-assemb-manifest_182'))
    meta_list = []
    row1 = ''
    id_check = []
    with open(os.path.join(ra_dir, 'f_all.csv')) as in_f:
        for line in in_f:
            if line.startswith('#SampleID'):
                row1 = line.rstrip().replace(',', '\t')
            else:
                lf = line.rstrip().split(',')
                lf[1] = str(secrets.choice(range(20, 70))) if lf[1] == 'NaN' else lf[1]
                tmp = '\t'.join(lf)
                if lf[0] in s_set:
                    id_check.append(lf[0])
                    meta_list.append(tmp)
    with open(os.path.join(ra_dir, 'ra-s-metadata_f.tsv'), 'w') as o_f:
        o_f.write(row1 + '\n')
        o_f.write('\n'.join(meta_list) + '\n')
        modify_man(os.path.join(ra_dir, 'ra-assemb-manifest_182'), id_check)


if __name__ == '__main__':
    main()
