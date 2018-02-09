#!/usr/bin/env python3.6
import os
import sys
import random
import argparse


def args_parse():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-t', '--type', type=str, choices=['p', 's'], default='s',
                        help='read type, p: paired-end, s: single-end, default is "s"')
    parser.add_argument('directory', help='directory of reads')
    parser.add_argument('-m', '--meta', action="store_true", help='To create pseudo-metadata when set "-m"')
    args = parser.parse_args()
    return args


def get_format(fq_list, abs, end_type):
    result_list = []
    #f_count = 0
    for fq in fq_list:
        td = fq[:fq.find('.')]
        if end_type == 's':
            #f_count += 1
            result_list.append('sample{},{}/{},forward'.format(td, abs, fq))
        else:
            if 'R1' in fq:
                #f_count += 1
                fq2 = fq.replace('R1', 'R2')
                result_list.append('sample{},{}/{},forward\nsample{},{}/{},reverse'.format(td, abs, fq, td, abs, fq2))
    return result_list


def pseudo_metadata(r_list, end_type):
    prefix = 'ra-s' if end_type == 's' else 'ra-p'
    response = ['high', 'mid', 'low']
    gender = ['female', 'male']
    drug = ['Etanercept', 'Adalimumab', 'Golimumab', 'Tocilizumab', 'Abatacept']
    with open('{}-metadata.tsv'.format(prefix), 'w') as out_f:
        out_f.write('#SampleID\tresponse\tgender\tdate\tdrug\n')
        for ele in r_list:
            s_id = ele.split(',')[0]
            res = random.choice(response)
            sex = random.choice(gender)
            day = random.randint(1, 301)
            trt = random.choice(drug)
            out_f.write('{}\t{}\t{}\t{}\t{}\n'.format(s_id, res, sex, day, trt))


def main(argv=None):
    if argv is None:
        argv = args_parse()
        print(argv)
        fq_list = os.listdir(argv.directory)
        fq_list.sort()
        abs_path = os.path.abspath(argv.directory)
        fpath, fname = os.path.split(abs_path)
        fbase, fext = os.path.splitext(fname)
        result_list = get_format(fq_list, abs_path, argv.type)
        if argv.meta:
            pseudo_metadata(result_list, argv.type)
        with open('ra-{}-manifest'.format(fbase), 'w') as out_f:
            out_f.write('sample-id,absolute-filepath,direction\n')
            out_f.write('\n'.join(result_list) + '\n')


if __name__ == '__main__':
    sys.exit(main())
