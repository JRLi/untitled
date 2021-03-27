#!/usr/bin/env python
import pandas as pd
import numpy as np
import os
dir_htseq = '/mount/ictr1/chenglab/jli/Projects/RNA_expression/htseq_out/hisat2/GSE133789'
dir_out = '/mount/ictr1/chenglab/jli/Projects/RNA_expression/TPM/drug_inhibit/GSE133789'
len_file = '/mount/ictr1/chenglab/jli/db/ensembl/saccharomyces/Saccharomyces_cerevisiae_len'
suf = '_Introns.counts.txt'


def mkdir(dir_in):
    if not os.path.isdir(dir_in):
        os.mkdir(dir_in)


def calculate_tpm_exon(df_in, df_g):
    df = df_in.copy()
    ss_e = df.exon.div(df_g.len_exon) * 1000
    ss_e.replace([np.inf, -np.inf], np.nan, inplace=True)
    e_tpm = ss_e.divide(ss_e.sum()) * 1000000
    return e_tpm


def calculate_tpm(df_in, df_g):
    df = df_in.copy()
    ss_e = df.exon.div(df_g.len_exon) * 1000
    ss_e.replace([np.inf, -np.inf], np.nan, inplace=True)
    ss_i = df.intron.div(df_g.len_intron) * 1000
    ss_i.replace([np.inf, -np.inf], 0, inplace=True)
    e_tpm = ss_e.divide(ss_e.sum() + ss_i.sum()) * 1000000
    i_tpm = ss_i.divide(ss_e.sum() + ss_i.sum()) * 1000000
    return e_tpm, i_tpm


def parse_count(file_list):
    df_len = pd.read_table(len_file, index_col=0)
    df_out_e = pd.DataFrame()
    df_out_ee = pd.DataFrame()
    df_out_i = pd.DataFrame()
    df_out_x = pd.DataFrame()
    df_count = pd.DataFrame()
    for f_i, f_e in zip(file_list[::2], file_list[1::2]):
        f_pre = f_i.strip().replace(suf, '')
        df_e = pd.read_table(os.path.join(dir_htseq, f_e), index_col=0, header=None, names=['exon'])
        df_i = pd.read_table(os.path.join(dir_htseq, f_i), index_col=0, header=None, names=['intron'])
        df_e = df_e.iloc[:-5, :]
        df_i = df_i.iloc[:-5, :]
        # df_i2 = df_i.fillna(0)
        dfs = pd.concat([df_e, df_i], 1)
        dfs['sum'] = dfs.sum(1)
        df_count[f_pre] = dfs.sum()
        dfs.fillna(0, inplace=True)
        tpm_e = calculate_tpm_exon(dfs, df_len)
        tpm_exon, tpm_intron = calculate_tpm(dfs, df_len)
        df_out_e[f_pre] = tpm_exon
        df_out_ee[f_pre] = tpm_e
        df_out_i[f_pre] = tpm_intron
        tpm_exon.index += '_exon'
        tpm_intron.index += '_intron'
        tpm_x = pd.concat([tpm_exon, tpm_intron])
        df_out_x[f_pre] = tpm_x
    df_out_e.dropna(how='all', inplace=True)
    df_out_ee.dropna(how='all', inplace=True)
    df_out_i.dropna(how='all', inplace=True)
    df_out_x.dropna(how='all', inplace=True)
    df_out_e = df_out_e.loc[(df_out_e != 0).any(axis=1)]
    df_out_ee = df_out_ee.loc[(df_out_ee != 0).any(axis=1)]
    df_out_i = df_out_i.loc[(df_out_i != 0).any(axis=1)]
    df_out_x = df_out_x.loc[(df_out_x != 0).any(axis=1)]
    df_out_e.to_csv(os.path.join(dir_out, 'TPM_exon'), sep='\t')
    df_out_ee.to_csv(os.path.join(dir_out, 'TPM_exon_nBYe'), sep='\t')
    df_out_i.to_csv(os.path.join(dir_out, 'TPM_intron'), sep='\t')
    df_out_x.to_csv(os.path.join(dir_out, 'TPM_all'), sep='\t')
    df_count.T.to_csv(os.path.join(dir_out, 'read_counts_summary'), sep='\t')


def main():
    f_list = os.listdir(dir_htseq)
    f_list.sort()
    mkdir(dir_out)
    parse_count(f_list)


if __name__ == '__main__':
    main()
