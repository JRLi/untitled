import os
import sys
import pandas as pd
import numpy as np
import scipy.stats as st
exo_f = 'D:/Project/lncRNA/exoRbase'


def open_df(in_path, direct='n'):
    fpath, fname = os.path.split(in_path)
    fbase, fext = os.path.splitext(fname)
    df = pd.read_csv(in_path, index_col=0) if fext == '.csv' else pd.read_table(in_path, index_col=0)
    if direct == 't':
        df = df.transpose()
    return df, fbase


def log2_transfer(df_in, p_c):
    dfi = df_in.copy()
    return np.log2(dfi + p_c)


def scipy_ttest_ind(s1, s2, var):
    return st.ttest_ind(s1, s2, equal_var=var)

def cn_t_test(dfn, dfd, var):
    dfo = pd.DataFrame(index=dfn.index)
    print(dfn.shape, dfd.shape)
    p_list = []
    for i in dfn.index:
        r = scipy_ttest_ind(dfn.ix[i, :], dfd.ix[i, :], var)
        p_list.append(r[1])
    dfo['P value'] = np.array(p_list)
    dfo['COAD_N'] = dfn.mean(1)
    dfo['COAD_T'] = dfd.mean(1)
    dfo['Fold-change T/N'] = dfd.mean(1) / dfn.mean(1)
    #dfo.to_csv(os.path.join(f_dir, dfb + '_t_test.csv'), na_rep=1)
    return dfo


dfc, dfb1 = open_df('D:/Project/lncRNA/tcga/COAD_Colon_adenocarcinoma_cancer_lc.csv')
dfd, dfc1 = open_df('D:/Project/lncRNA/tcga/COAD_Colon_adenocarcinoma_adjacent_normal_lc.csv')
r = scipy_ttest_ind(dfc.ix[0, :], dfd.ix[0, :], False)

f_dir = 'D:/Project/lncRNA/data'
f_s = 'Samples_combined_lncRNA_TPM.txt'
n_p = 'N'
d_p = 'CRC'
#  CRNDE/FENDRR and CCAT1/FENDRR
g_list = ['CCAT1', 'CRNDE', 'FENDRR']
df1, dfb = open_df(os.path.join(f_dir, f_s))
df1 = df1.ix[g_list]
n_c = [col for col in df1 if col.startswith(n_p)]
d_c = [col for col in df1 if col.startswith(d_p)]
dfn = df1[n_c]
dfc = df1[d_c]
dfo = pd.concat([dfn, dfc], 1)
#dfo.to_csv(os.path.join(f_dir, 'exoRbase_NC_3gene.csv'))
dfg = dfo + 0.1
s1 = dfg.ix[0]/dfg.ix[2]
s1.name = 'CCAT1/FENDRR_pc0.1'
s2 = dfg.ix[1]/dfg.ix[2]
s2.name = 'CRNDE/FENDRR_pc0.1'
dfo = dfo.append(s1)
dfo = dfo.append(s2)
#print(dfo)
#dfo.T.to_csv(os.path.join(f_dir, 'exoRbase_NC_3gene_pc0.1T.csv'))
dfon = dfo[n_c]
dfoc = dfo[d_c]
df_f = cn_t_test(dfon, dfoc, False)
print(df_f)
#df_f.to_csv(os.path.join(f_dir, 'coad_3gene_ttest.csv'))


def g21_e2s(path):
    with open(path) as in_f:
        e2s_dict, s2e_dict = {}, {}
        for line in in_f:
            lf = line[1:].rstrip().split('|')
            e2s_dict[lf[0]] = lf[5]
            s2e_dict[lf[5]] = lf[0]
        print(len(s2e_dict), len(e2s_dict))
        print(e2s_dict.get('-'))
        print(s2e_dict.get('-'))
        return e2s_dict


def exoRbase_it_COAD():
    with open('D:/Project/lncRNA/tcga/GSE_COAD_cvn_lc_ttest.csv') as in_fg, \
            open('D:/Project/lncRNA/tcga/Samples_combined_lncRNA_TPM_t_test.csv') as in_fe, \
            open('D:/Project/lncRNA/tcga/TCGA_COAD_Colon_adenocarcinoma_cvn_lc_t_test.csv') as in_ft, \
            open('D:/Project/lncRNA/tcga/GSE_COAD_exoRbase.csv', 'w') as out_feg, \
            open('D:/Project/lncRNA/tcga/TCGA_COAD_exoRbase.csv', 'w') as out_fet:
        next(in_fe)
        next(in_fg)
        next(in_ft)
        out_feg.write('lncRNA,P_value_gse,N_mean,T_mean,Fold-change T/N,P_value_exo,N_mean,T_mean,Fold-change T/N\n')
        out_fet.write('lncRNA,P_value_tcga,N_mean,T_mean,Fold-change T/N,P_value_exo,N_mean,T_mean,Fold-change T/N\n')
        ed = dict([x.rstrip().split(',', 1) for x in in_fe])
        for line in in_fg:
            lf = line.rstrip().split(',', 1)
            if lf[0] in ed:
                out_feg.write(','.join(lf) + ',' + ed.get(lf[0]) + '\n')
        for line in in_ft:
            lf = line.rstrip().split(',', 1)
            if lf[0] in ed:
                out_fet.write(','.join(lf) + ',' + ed.get(lf[0]) + '\n')


def exoRbase_it_LIHC():
    with open('D:/Project/lncRNA/tcga/GSE63420_liver_cvn_lc_t_test.csv') as in_fg, \
            open('D:/Project/lncRNA/tcga/Samples_combined_lncRNA_TPM_t_test.csv') as in_fe, \
            open('D:/Project/lncRNA/tcga/TCGA_LIHC_Liver_hepatocellular_carcinoma_cvn_lc_t_test.csv') as in_ft, \
            open('D:/Project/lncRNA/tcga/GSE63420_LIHC_exoRbase.csv', 'w') as out_feg, \
            open('D:/Project/lncRNA/tcga/TCGA_LIHC_exoRbase.csv', 'w') as out_fet:
        next(in_fe)
        next(in_fg)
        next(in_ft)
        out_feg.write('lncRNA,P_value_gse,N_mean,T_mean,Fold-change T/N,P_value_exo,N_mean,T_mean,Fold-change T/N\n')
        out_fet.write('lncRNA,P_value_tcga,N_mean,T_mean,Fold-change T/N,P_value_exo,N_mean,T_mean,Fold-change T/N\n')
        ed = dict([x.rstrip().split(',', 1) for x in in_fe])
        for line in in_fg:
            lf = line.rstrip().split(',', 1)
            if lf[0] in ed:
                out_feg.write(','.join(lf) + ',' + ed.get(lf[0]) + '\n')
        for line in in_ft:
            lf = line.rstrip().split(',', 1)
            if lf[0] in ed:
                out_fet.write(','.join(lf) + ',' + ed.get(lf[0]) + '\n')


def e2c_lihc(path):
    g_l = ['AFP', 'AHSG', 'ALB', 'APOH', 'FABP1', 'FGB', 'FGG', 'GPC3', 'RBP4', 'TF']
    df1, dfb = open_df(path)
    print(df1.shape)
    df2 = df1.ix[g_l, :]
    print(df2.shape)
    df2.to_csv(os.path.join(exo_f, 'mRNA_TPM_10genes_NH.csv'))


def main():
    #e2s_dict = g21_e2s('D:/Project/lncRNA/tcga/g21_ln_list')
    #exoRbase_it_COAD()
    #exoRbase_it_LIHC()
    e2c_lihc(os.path.join(exo_f, 'Samples_combined_mRNA_TPM.txt'))


if __name__ == '__main__':
    main()