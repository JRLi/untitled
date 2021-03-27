#!/usr/bin/env python
import pandas as pd
import os
gtf_exo = '/home/233549/jli/db/ensembl/saccharomyces/Saccharomyces_cerevisiae.R64.consExons.gtf'
gtf_int = '/home/233549/jli/db/ensembl/saccharomyces/Saccharomyces_cerevisiae.R64.Introns.gtf'
#gtf_exo = '/home/233549/jli/db/ensembl/gtf/Homo_sapiens.GRCh38.99.chr.consExons.gtf'
#gtf_int = '/home/233549/jli/db/ensembl/gtf/Homo_sapiens.GRCh38.99.chr.Introns.gtf'
# gtf_exo = 'C:/Users/fanic/Downloads/tmp/Homo_sapiens.GRCh37.87.chr.consExons.gtf'
# gtf_int = 'C:/Users/fanic/Downloads/tmp/Homo_sapiens.GRCh37.87.chr.Introns.gtf'
out_len = 'Saccharomyces_cerevisiae_len'


def df_open(path_in):
    p_d, f_n = os.path.split(path_in)
    f_p, f_s = os.path.splitext(f_n)
    df = pd.read_csv(path_in, header=None) if f_s == '.csv' else pd.read_table(path_in, header=None)
    return df, p_d, f_p


def df_modify(df_in):
    df = df_in.copy()
    df[8] = df[8].str.replace('gene_id "', '')
    df[8] = df[8].str.replace('";', '')
    df[9] = df[4] - df[3] + 1
    df.sort_values(8, inplace=True)
    df.drop([1, 5, 7], axis=1, inplace=True)
    return df


def main():
    dfe, dfe_p, dfe_f = df_open(gtf_exo)
    dfi, dfi_p, dfi_f = df_open(gtf_int)
    dfe = df_modify(dfe)
    dfi = df_modify(dfi)
    # dfe.to_csv(os.path.join(dfe_p, '{}_withLength'.format(dfe_f)), sep='\t', header=False, index=False)
    # dfi.to_csv(os.path.join(dfi_p, '{}_withLength'.format(dfi_f)), sep='\t', header=False, index=False)
    dfe = dfe.groupby(8)[[9]].sum()
    dfi = dfi.groupby(8)[[9]].sum()
    # dfe.to_csv(os.path.join(dfe_p, '{}_G2length'.format(dfe_f)), sep='\t', header=False)
    # dfi.to_csv(os.path.join(dfi_p, '{}_G2length'.format(dfi_f)), sep='\t', header=False)
    dfe.columns = ['len_exon']
    dfi.columns = ['len_intron']
    dfs = pd.concat([dfe, dfi], 1).fillna(0)
    dfs.to_csv(os.path.join(dfe_p, out_len), sep='\t')


if __name__ == '__main__':
    main()
