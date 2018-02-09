import xmltodict
import os
import pandas as pd
import numpy as np
from sklearn.preprocessing import Imputer
r_p = 'E:/StringTemp/Project_Rice/'


def open_df(in_path, direct='n'):
    fpath, fname = os.path.split(in_path)
    fbase, fext = os.path.splitext(fname)
    df = pd.read_csv(in_path, index_col=0) if fext == '.csv' else pd.read_table(in_path, index_col=0)
    print('{}: Transpose: {}'.format(fbase, direct))
    if direct == 't':
        df = df.transpose()
    return df, fbase


def sk_imputation(df_in):
    df = df_in.copy()
    imr = Imputer(missing_values='NaN', strategy='mean', axis=0)
    imr = imr.fit(df)
    arr = imr.transform(df)
    return arr


def filter_novel(n_path, f_path, threshold=2):
    dfn, dfn_n = open_df(n_path)
    dff, dff_n = open_df(f_path)
    idx = dff[dff.ix[:, 0] >= threshold].index
    dfn = dfn.ix[idx]
    return dfn


def rpm(df_input):
    df = df_input.copy()
    return (df / df.sum()) * 1000000


def quantileNormalize(df_input):
    df = df_input.copy()
    #compute rank
    dic = {}
    for col in df:
        dic.update({col : sorted(df[col])})
    sorted_df = pd.DataFrame(dic)
    rank = sorted_df.mean(axis = 1).tolist()
    #sort
    for col in df:
        t = np.searchsorted(np.sort(df[col]), df[col])
        df[col] = [rank[i] for i in t]
    return df


df1, df_b = open_df(os.path.join(r_p, 'nm_df129_q.csv'))
arr1 = sk_imputation(df1)
df = pd.DataFrame(arr1, index=df1.index, columns=df1.columns)
df.to_csv(os.path.join(r_p, 'nm_df129_q_i.csv'))


"""
dfn = filter_novel(os.path.join(r_p, 'novel_c2.csv'), os.path.join(r_p, 'novel_c2_n0.csv'), 10)
dfb, dfb_n = open_df(os.path.join(r_p, 'blast_df_0.8_0.8_0.9_nm.csv'))
dff = pd.concat([dfb, dfn])
#dff.to_csv(os.path.join(r_p, 'rice_nov_mir_wm.csv'))
#dff = quantileNormalize(dff)
dff = rpm(dff)
#dff.to_csv(os.path.join(r_p, 'quant_nm.csv'))
dfp, dfp_n = open_df(os.path.join(r_p, 'phenotype.csv'), 't')
dfp = dfp.add_suffix('p')
idx = dfp.columns.intersection(dff.columns)
df = pd.concat([dff, dfp], join='inner')
#df.T.to_csv(os.path.join(r_p, 'nm_df129_r.csv'), na_rep='NA')
"""

'''
d_r = 'D:/Project/circ_miRNA/search_result'

ld = os.listdir(d_r)
with open('D:/Project/circ_miRNA/clinical_p', 'w', encoding = 'utf8') as out_f:
    for f in ld:
        with open(os.path.join(d_r, f), encoding = 'utf8') as in_f:
            for line in in_f:
                if 'mir-' in line or 'miR-' in line:
                    print(line)
                    out_f.write('{}\t{}\n'.format(f, line.rstrip()))
        #od = xmltodict.parse(in_f.read())
        #print(od['clinical_study'].keys())
        #for k, v in od['clinical_study'].items():
            #print(k, v, sep='\n')
'''