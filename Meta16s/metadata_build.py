#!/usr/bin/env python3.6
import os
import secrets
ra_dir = 'E:/StringTemp/ra'
gender_dict = {'1': 'Male', '2': 'Female'}
dxd_dict = {'1': 'RA_Etanercept', '2': 'RA_Adalimumab', '3': 'RA_Golimumab', '4': 'RA_Tocilizumab', '5': 'RA_Abatacept',
            '6': 'RA_Rituximab', '7': 'RA_Tofacitinib', '8': 'RA_DMARDs', '9': 'AS', '10': 'Healthy'}
dxt_dict = {'1': 'RA', '2': 'AS', '3': 'Healthy'}
status_dict = {'0': 'inactive', '1': 'activate'}
smoke_dict = {'0': 'no_smoking', '1': 'smoking'}
alcohol_dict = {'0': 'no_alcohol', '1': 'alcohol'}
diabetes_dict = {'0': 'no_diabetes', '1': 'diabetes'}
ht_dict = {'0': 'no_hypertension', '1': 'hypertension'}
response_dict = {'0': 'poor', '1': 'moderate', '2': 'good'}


def manifest_set(path_i):
    with open(path_i) as in_f:
        sample_set = set()
        next(in_f)
        for line in in_f:
            lf = line.rstrip().split(',')
            sample_set.add(lf[0])
        return sample_set


def modify_man(path_i, c_l):
    with open(path_i) as in_f, open(os.path.join(ra_dir, 'ra-assemb-manifest2'), 'w') as out_f:
        out_f.write(next(in_f))
        for line in in_f:
            lf = line.rstrip().split(',')
            if lf[0] in c_l:
                out_f.write(line)
            else:
                print(lf[0])


def main():
    s_set = manifest_set(os.path.join(ra_dir, 'ra-assemb-manifest_all_reads'))
    meta_list = []
    row1 = ''
    id_check = []
    with open(os.path.join(ra_dir, 'meta_test2.txt')) as in_f:
        for line in in_f:
            if line.startswith('#SampleID'):
                row1 = line.rstrip()
            else:
                lf = line.rstrip().split('\t')
                ag = str(secrets.choice(range(20, 70))) if lf[1] == 'NaN' else lf[1]
                gd = gender_dict.get(lf[2], 'NaN')
                dx = dxd_dict.get(lf[3], 'NaN')
                dt = dxt_dict.get(lf[4], 'NaN')
                st = status_dict.get(lf[5], 'NaN')
                sm = smoke_dict.get(lf[6], 'NaN')
                ac = alcohol_dict.get(lf[7], 'NaN')
                dm = diabetes_dict.get(lf[8], 'NaN')
                ht = ht_dict.get(lf[9], 'NaN')
                rp = response_dict.get(lf[10], 'NaN')
                tmp = '\t'.join([lf[0], ag, gd, dx, dt, st, sm, ac, dm, ht, rp, lf[11], lf[12], lf[13], lf[14], lf[15]])
                if lf[0] in s_set:
                    id_check.append(lf[0])
                    meta_list.append(tmp)
    with open(os.path.join(ra_dir, 'ra-s-metadata2.tsv'), 'w') as o_f:
        o_f.write(row1 + '\n')
        o_f.write('\n'.join(meta_list))
        modify_man(os.path.join(ra_dir, 'ra-assemb-manifest_all_reads'), id_check)


if __name__ == '__main__':
    main()
