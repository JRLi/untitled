#!/usr/bin/env python
import os
from collections import Counter
path = 'D:/Project/circ_miRNA'


def get_c_tcga(file):
    fbase, fext = os.path.splitext(file)
    with open(os.path.join(path, file)) as in_f:
        next(in_f)
        c_list = [line.split(',')[2].split('-')[3][0:2] for line in in_f]
        s_counter = Counter(c_list)
        return fbase, len(c_list), s_counter


def get_c_target(file):
    fbase, fext = os.path.splitext(file)
    with open(os.path.join(path, file)) as in_f:
        next(in_f)
        s_list = []
        for line in in_f:
            case = line.split(',')[2].split('-')
            if len(case) < 4:
                print(fbase, line)
            else:
                check = case[3][0:2]
                s_list.append(check)
        s_counter = Counter(s_list)
        return fbase, len(s_list), s_counter


def main():
    file_list = os.listdir(path)
    with open(path + '/check_gdc', 'w') as gdc_o:
        gdc_o.write('project\tall_case\ttumor\tcontrol\n')
        for file in file_list:
            if file.startswith('TCGA'):
                proj, all_s, sample_counter = get_c_tcga(file)
            elif file.startswith('TARGET'):
                proj, all_s, sample_counter = get_c_target(file)
            else:
                continue
            print('{}\t{}'.format(proj, all_s), end='\t')
            gdc_o.write('{}\t{}\t'.format(proj, all_s))
            t_count, c_count, c_list = 0, 0, []
            for k, v in sample_counter.items():
                c_list.append('{}: {}'.format(k, v))
                if int(k) <= 9:
                    t_count += v
                else:
                    c_count += v
            gdc_o.write('{}\t{}\n'.format(t_count, c_count))
            print('\t'.join(c_list))

if __name__ == '__main__':
    main()
