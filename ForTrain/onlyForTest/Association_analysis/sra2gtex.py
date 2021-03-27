#!/usr/bin/env python
import pandas as pd
import os
dir_exp = '/mount/ictr1/chenglab/jli/Projects/RNA_expression/TPM/new'
emt = '/mount/ictr1/chenglab/jli/Projects/RNA_expression/ID_annotation/E-MTAB-5214-sdrf'


def srr_gtex_dict(in_file):
    s2g_dict, g2o_dict = dict(), dict()
    with open(in_file, 'r') as in_f:
        next(in_f)
        for line in in_f:
            line_ls = line.rstrip().split('\t')
            s2g_dict[line_ls[14]] = line_ls[0]
            g2o_dict[line_ls[0]] = line_ls[5]
        return s2g_dict, g2o_dict


def df_mean_columns(df_input):
    gp = df_input.groupby(df_input.columns, axis=1)
    return gp.mean()


def df_modify(files, s2g, g2o, col_s):
    for f in files:
        if f in ['TPM_intron_pro_rc', 'TPM_intron_pro_tpm', 'intron_pro_rc', 'intron_pro_tpm']:
        # if f in ['stab_99', 'stab_01']:
            print(f)
            df = pd.read_table(os.path.join(dir_exp, f), index_col=0)
            if f.startswith('stab'):
                col = [x.split('.')[0] for x in df.columns]
                df.columns = col
            df = df[col_s]
            df = df.loc[(df != 0).any(axis=1), (df != 0).any(axis=0)]  # drop with all zero
            cd_list = [col for col in df.columns if s2g.get(col) is None]
            df.drop(cd_list, 1, inplace=True)
            df.rename(s2g, axis=1, inplace=True)
            df.sort_index(axis=1, inplace=True)
            print('before mean: ', len(df.columns))
            with open(os.path.join(dir_exp, 'g2o'), 'w') as out_f:
                for g in df.columns:
                    out_f.write('{}\t{}\n'.format(g, g2o.get(g)))
            df2 = df.copy()
            df2.columns = ['-'.join(x.split('-')[:2]) for x in df2.columns]
            df2 = df_mean_columns(df2)
            print('after mean: ', len(df2.columns))
            #df.to_csv(os.path.join(dir_exp, '{}_gIDa'.format(f)), sep='\t')
            df2.to_csv(os.path.join(dir_exp, '{}_gIDd_mean'.format(f)), sep='\t')
            # df.to_csv(os.path.join(dir_exp, 'TPM_gtexID_full_{}'.format(f.split('_')[1])), sep='\t')
            # df2.to_csv(os.path.join(dir_exp, 'TPM_gtexID_short_mean_{}'.format(f.split('_')[1])), sep='\t')
            # df.to_csv(os.path.join(dir_exp, 'mx_{}_gIDa'.format(f.split('.')[3])), sep='\t')
            # df2.to_csv(os.path.join(dir_exp, 'mx_{}_gIDd_mean'.format(f.split('.')[3])), sep='\t')


def main():
    sra2gext, gext2organism = srr_gtex_dict(emt)
    f_list = os.listdir(dir_exp)
    df = pd.read_table(os.path.join(dir_exp, 'read_counts_summary'), index_col=0)
    print(df.shape)
    id1 = df[df.exon < 3000000].index
    print('less than 3M align')
    print(len(id1))
    id2 = [x for x in df.index if x.startswith('SRR21')]
    print('remove abnormal {}'.format(len(id2)))
    idx = id1.union(id2)
    df.drop(idx, inplace=True)
    df_modify(f_list, sra2gext, gext2organism, df.index)


if __name__ == '__main__':
    main()
