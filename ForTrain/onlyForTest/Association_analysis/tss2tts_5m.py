#!/usr/bin/env python
import os, sys
import pandas as pd
import numpy as np
from _collections import defaultdict
ensGene = '/mount/amos1/home/jli/Project_tmp/phg000520.v2.GTEx_MidPoint_Imputation.genotype-calls-vcf.c1/ensGene'
dosage_f = '/mount/amos1/home/jli/Project_tmp/phg000520.v2.GTEx_MidPoint_Imputation.genotype-calls-vcf.c1/' \
           'dosage-out1-rm-iloc-inter-df'
out_dir = '/mount/amos1/home/jli/Project_tmp/phg000520.v2.GTEx_MidPoint_Imputation.genotype-calls-vcf.c1/trans2SNV_out'


def df_open(path_in, i_col, dt):
    p_d, f_n = os.path.split(path_in)
    f_p, f_s = os.path.splitext(f_n)
    df = pd.read_csv(path_in, index_col=i_col, dtype=dt) if f_s == '.csv' \
        else pd.read_table(path_in, index_col=i_col, dtype=dt)
    return df, p_d, f_p


def range_5m(df_in, rg):
    df = df_in.copy()
    df.txStart = df.txStart.apply(lambda x: x - rg if (x - rg > 0) else 0)
    df.txEnd = df.txEnd + rg
    df.chrom = df.chrom.replace(regex='chr', value='')
    return df


def main():
    rg = sys.argv[1]
    if rg is None:
        sys.exit()
    d_type = {'chrom': str, 'txStart': int, 'txEnd': int}
    df1, d1_p, d1_b = df_open(ensGene, 1, d_type)
    df2, d2_p, d2_b = df_open(dosage_f, 0, d_type)
    dt_dict = defaultdict(dict)
    for x in df2.index:
        dt_dict[x.split('_')[0]][x] = int(x.split('_')[1])
    # snv_dict = {k: k.split('_')[:2] for k in df2.index}
    df1 = df1[df1.columns[1:5]]
    df1 = range_5m(df1, float(rg))
    chr_set = set(df1.chrom)
    for chr_n in chr_set:
        chr_dict = dt_dict.get(chr_n)
        if chr_dict is not None:
            print(chr_n)
            target_snv = []
            df2 = df1.loc[df1.chrom == chr_n]
            for i in range(len(df2.index)):
                snv_list = []
                ss = df2.loc[df2.index[i]]
                for k, v in chr_dict.items():
                    if ss.txStart <= v <= ss.txEnd:
                        snv_list.append(k)
                target_snv.append(' '.join(snv_list))
            # df2['cisSNV'] = np.array(target_snv)
            df2['SNV'] = target_snv
            df2.to_csv(os.path.join(out_dir, 'tanscript2SNV_chr{}_{}'.format(chr_n, rg)), sep='\t')


if __name__ == '__main__':
    main()
