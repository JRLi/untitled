#!/usr/bin/env python
import scipy.stats.mstats as mt
import pandas as pd
import os
lincs_exp_path = './'


def get_probe2gene_dict(probe2geneID_path):
    probe2gene_dict = {}  # value is list type
    with open(probe2geneID_path, 'r') as p2g_file:
        for line in p2g_file:
            lf = line.replace('\r\n', '').replace('\n', '').replace('\"', '').split('\t')
            g_list = probe2gene_dict.get(lf[0], [])
            g_list.append(lf[1])
            probe2gene_dict[lf[0]] = g_list
        i, j = 0, 0
        for k, v in probe2gene_dict.items():
            if len(v) > 1:
                print(k, v)
                i += 1
            else:
                j += 1
        print(i, j)
        print(len(probe2gene_dict))
        return probe2gene_dict


def out_cmap_with_geneid(cmap_path, probe2gene_dict):
    with open(cmap_path, 'r') as cmap_in, open('./rankMatrix_covertGID.txt', 'w') as cmap_out, \
            open('./rankMatrix_noGID.txt', 'w') as cmap_no, open('./cmap_wP2G.txt', 'w') as id_doc:
        for line in cmap_in:
            if line.startswith('probe_id'):
                cmap_out.write(line.replace('\r\n', '').replace('\n', '').replace('\"', '') + '\n')
                cmap_no.write(line.replace('\r\n', '').replace('\n', '').replace('\"', '') + '\n')
            else:
                lf = line.replace('\r\n', '').replace('\n', '').replace('\"', '').split('\t', maxsplit=1)
                g_id_list = probe2gene_dict.get(lf[0])
                if g_id_list is None:
                    cmap_no.write(line.replace('\r\n', '').replace('\n', '').replace('\"', '') + '\n')
                else:
                    for g_id in g_id_list:
                        cmap_out.write('{}\t{}\n'.format(g_id, lf[1]))
                        id_doc.write('{}\t{}\n'.format(lf[0], g_id))


def process_gdsc_exp():
    ensembl2gene_dict = {}
    gdsc_exp_path = 'D://Project/drs/gdsc/sanger1018_brainarray_ensemblgene_rma.txt'
    hgnc_path = 'D://Project/drs/gdsc/cmap/hgnc_4_170407'
    cmap_dir = 'D://Project/drs/gdsc/cmap/'
    a, b, c = 0, 0, 0
    with open(hgnc_path, 'r') as hgnc_in:
        for line in hgnc_in:
            if line.startswith('hgnc_id'):
                continue
            lf = line.replace('\r\n', '').replace('\n', '').replace('\"', '').split('\t')
            if lf[3] != '':
                a += 1
                if lf[2] != '':
                    b += 1
                    if lf[3] in ensembl2gene_dict.keys():
                        print(lf[3], lf[2], ensembl2gene_dict.get(lf[3]))
                    ensembl2gene_dict[lf[3]] = lf[2]
                else:
                    c += 1
        print(len(ensembl2gene_dict))
        print(a, b, c)
    with open (gdsc_exp_path, 'r') as gdsc_in, open(cmap_dir + 'gdsc_w_exp.txt', 'w') as out1, \
            open(cmap_dir + 'gdsc_no_exp.txt', 'w') as out2:
        for line in gdsc_in:
            if line.startswith('ensembl_gene'):
                out1.write(line.replace('\r\n', '').replace('\n', '').replace('\"', '') + '\n')
                out2.write(line.replace('\r\n', '').replace('\n', '').replace('\"', '') + '\n')
            else:
                lf = line.replace('\r\n', '').replace('\n', '').replace('\"', '').split('\t', maxsplit=1)
                gid = ensembl2gene_dict.get(lf[0])
                if gid is None:
                    out2.write(line.replace('\r\n', '').replace('\n', '').replace('\"', '') + '\n')
                else:
                    out1.write('{}\t{}\n'.format(gid, lf[1]))


def openDF(in_path):
    fpath, fname = os.path.split(in_path)
    fbase, fext = os.path.splitext(fname)
    df = pd.read_csv(in_path, index_col=0) if fext == '.csv' else pd.read_table(in_path, index_col=0)
    return df, fbase


def df_mean_index(df_input):
    gp = df_input.groupby(df_input.index)
    return gp.mean()


def z_transfer_mode_select(df_input, mode='n'):
    df = df_input.copy()
    # mean and std track column from row index 0 to end, so axis use default (index 0).
    avg = df.mean()
    std = df.std()
    for k in range(len(df.columns)):
        df.iloc[:, k] = (df.iloc[:, k] - avg[k]) / std[k]
    if mode == 'n':
        return df
    else:
        return -df


def scipy_z_transfer(df_input):
    df = df_input.copy()
    for k in range(len(df.columns)):
        df.iloc[:, k] = mt.zscore(df.iloc[:, k], ddof=1)
    return df


def main():
    # process_gdsc_exp()
    # p2id_dict = get_probe2gene_dict('./LINCS_GPL20573_probe2geneid.txt')
    # out_cmap_with_geneid('./rankMatrix.txt', p2id_dict)
    df1, file_base = openDF('./rankMatrix_covertGID.txt')
    print(df1.shape)
    df1 = df_mean_index(df1)
    print(df1.shape)
    df1 = z_transfer_mode_select(df1, 'r')  # reverse Z
    print(df1.shape)
    df1.to_csv('./CMAP_rank_reverseZscore.csv')

if __name__ == '__main__':
    main()