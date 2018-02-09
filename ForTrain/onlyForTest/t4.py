import os
import pandas as pd
import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt

def reverse_complement(line, seq = False):
    complement_dict = {'a':'t','t':'a','c':'g','g':'c','A':'T','T':'A','C':'G','G':'C'}
    r_line = line[::-1]
    if seq:
        rc_line = "".join(complement_dict.get(bp, bp) for bp in r_line)
        return rc_line
    else:
        return r_line


def insert_end(line, index):
    return "".join([line[j:j+index]+"\n" for j in range(0, len(line), index)])

read2 = 'AGGTTCCGGATAAGTAAGAGCCT'
seqr = reverse_complement(read2, True)
print('test', seqr)


def open_df(in_path, direct='n'):
    fpath, fname = os.path.split(in_path)
    fbase, fext = os.path.splitext(fname)
    df = pd.read_csv(in_path, index_col=0) if fext == '.csv' else pd.read_table(in_path, index_col=0)
    print('{}: Transpose: {}'.format(fbase, direct))
    if direct == 't':
        df = df.transpose()
    return df, fbase


def scipy_corr(s1, s2, corr_mode):
    ixs = s1.index.intersection(s2.index)
    if corr_mode == 'pearson':
        return st.pearsonr(s1[ixs], s2[ixs])
    elif corr_mode == 'kendall':
        return st.kendalltau(s1[ixs], s2[ixs])
    elif corr_mode == 'spearman':
        return st.spearmanr(s1[ixs], s2[ixs])


df1, df1_base = open_df('D:/Project/circ_miRNA/forLi_mRNA_andMiRNA/test.txt')
print(df1)
print(df1.ix['TCGA-BLCA'])
print(scipy_corr(df1.ix['TCGA-BLCA'], df1.ix['TCGA-BLCA'], 'pearson'))

df2, df2_base = open_df('D:/Project/circ_miRNA/forLi_mRNA_andMiRNA/test_mi.csv')
print(df2)
print(df2.columns)
c_list = ['-'.join(x.split('-')[0:4]) for x in df2.columns]
print(c_list)
df2.columns = ['-'.join(x.split('-')[0:4]) for x in df2.columns]
print(df2)
ii1 = np.where(df2.iloc[0, :] < 10000)
ii2 = np.where(df2.iloc[4, :] < 10000)
print(ii1)
print(ii2)
print(type(df2.ix['hsa-let-7c', ii1[0]]))
s1 = df2.ix['hsa-let-7f-1', ii1[0]]
s2 = df2.ix['hsa-let-7f-2', ii2[0]]
print(s1.index)
print(s2.index)
ixs = s1.index.intersection(s2.index)
print(ixs)
s1 = s1[ixs]
s2 = s2[ixs]
print(s2)
#print(df2.ix['hsa-let-7b', [0, 3, 4]])

dict1 = {'g1':{'1': 1, '2': 2, '3':3}, 'g2':{'1':2, '3':'null', '4': 5}, 'g3':{'0': 0, '5': 3}, 'g4':{'A': 2}}
df3 = pd.DataFrame(dict1)
print(df3)

df4, df_b = open_df('D:/Project/circ_miRNA/forLi_mRNA_andMiRNA/test2.csv')
print(df4)
print(df4.ix[1028])


a = pd.Series([1,2,3], index=[1,1,2])
print(a)
print(a.index.duplicated())
a = a[~a.index.duplicated(keep='first')]
#a = a[a.index.drop_duplicates()]
print(a)

s1 = pd.Series([1, 2, 3, 4, 5, 6, 7, 8, 9], index = [1, 2, 3, 4, 5, 6, 7, 8, 9])
s2 = pd.Series([1, 3, 2, 4, 6, 5, 7, 9], index = [1, 3, 2, 4, 6, 5, 7, 9])
s3 = pd.Series([10, 30, 20, 40, 60, 50, 70, 90], index = [1, 3, 2, 4, 6, 5, 7, 9])
s4 = pd.Series([10, 30, 20, 40, 60, 50, 70, 90], index = [1, 2, 3, 4, 5, 6, 7, 8])
print(scipy_corr(s1, s2, 'pearson'))
print(scipy_corr(s1, s3, 'pearson'))
print(scipy_corr(s1, s4, 'pearson'))

dict1 = {'g1':{'1': 1, '2': 2, '3':3}, 'g2':{'1':2, '4': 5}, 'g3':{'0': 0, '5': 3}, 'g4':{'A': 2}}
df3 = pd.DataFrame(dict1)
df3 = df3.fillna(0)
print(df3)
#print(df3.ix[1])
df4 = df3[df3.ix[:,1] > 0]
#df5 = df3.ix[df3.ix[1, :] > 0]
#iim = np.where(df3.iloc[1, :] == 0)
#print(iim)
print(df4)
#print(df5)

df1, db = open_df('D:/Project/circ_miRNA/forLi_mRNA_andMiRNA/BRCA_hsa-mir-9-2_g1028_all.csv')
print(df1.iloc[0:5, :])
print(df1.corr())

#df1.plot.scatter(x='g1028', y='hsa-mir-9-2')
#plt.show()
#plt.matshow(df1.corr())
#plt.show()
import seaborn as sns; sns.set(color_codes=True)
ax = sns.regplot(x='g1028', y='hsa-mir-9-2', data=df1)
#plt.show()

list1 = [2, 12, 8, 0, 9]
list1.sort()
print(list1)
print(list1[0], list1[-1])

from sklearn import datasets
iris = datasets.load_iris()
df = pd.DataFrame(iris.data, columns=iris.feature_names)
print(df.head())

df1, b1 = open_df('D:/Project/orchid/all_new/test/t1.csv')
df2, b2 = open_df('D:/Project/orchid/all_new/test/t2.csv')
df3, b3 = open_df('D:/Project/orchid/all_new/test/t3.csv')
print(df1.ix[:, :1])    # get DF
print(df1.ix[:, 0])     # get Series

rdf = pd.concat([df1.ix[:, :1], df2.ix[:, :1], df3.ix[:, :1]], axis=1)
rdf2 = df1.ix[:, :1] + df2.ix[:, :1] + df3.ix[:, :1]

print(rdf)
print(rdf2)     # !!!!!!!!!not same with concat
rdf = rdf.fillna(0)
print(rdf)


import secrets
list1 = ['aaa', 'bbb', 'ccc']
print(secrets.choice(list1))

print(df1.T)
print(df2.T)
df3 = pd.concat([df1.T, df2.T], join='inner')
print(df3)