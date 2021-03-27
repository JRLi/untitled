#! /usr/bin/env python
import os, sys
from collections import defaultdict, OrderedDict
ensGene = 'C:/Users/fanic/Downloads/tmp/ensGene'
out_dir = 'C:/Users/fanic/Downloads/tmp/ens_t2g'


def prepare_output_dir(output_dir):
    if os.path.exists(output_dir):
        pass
    else:
        os.mkdir(output_dir)


def f2dict(i_path):
    with open(i_path, 'r') as in_f:
        chr_t2g = defaultdict(lambda: OrderedDict())
        next(in_f)
        for line in in_f:
            lss = line.rstrip().split('\t')
            chr_n = lss[2].replace('chr', '')
            if len(chr_n) <= 2:
                chr_t2g[chr_n][lss[1]] = lss[12]
        return chr_t2g


def out_dict(i_dict):
    for chr_n, t2g in i_dict.items():
        with open(os.path.join(out_dir, 't2g_chr{}'.format(chr_n)), 'w') as out_f:
            for t, g in i_dict.get(chr_n).items():
                out_f.write('{}\t{}\n'.format(t, g))


def main():
    prepare_output_dir(out_dir)
    chr2t2g = f2dict(ensGene)
    out_dict(chr2t2g)


if __name__ == '__main__':
    sys.exit(main())
