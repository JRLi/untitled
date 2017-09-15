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

seq_pp2 = 'CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT'
rev_pp2 = 'AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG'   # SRR3234010_1.fastq have 10698296
seq_pp1 = 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'
rev_pp1 = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'  # some R2 file contain this, like SRR3233954_2.fastq
seq_ill_uni = 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'
seq_3 = 'CGCTCTTGGGAACACGT' # SRR3234010_1.fastq have 9711954 at 4st nt
seq_2_pp2 = 'AGAGGACCTGAACAACGTGTTCCCAAGA'
read1='AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'
read2 = 'CGCCTTGGCCGTACAGCAGCCTCTTA'
seqr = reverse_complement(read2, True)
print('test', seqr)

barcode1 = 'TGACTA'
barcode2 = reverse_complement(barcode1, True)
readR1 = 'NAATGACTAGCAGGGTCAGGGTTCTGGATATCACACAGAGGTAGGTGGCTGAATCACTGA'
readR2 = 'AGCAGGAACAAGACTATTTGTTAAAGCAAATATCCAGAANCCTGACCCTGCTAGTCACTC'
print(barcode2)
print(readR1.find(barcode1))
print(readR2.rfind(barcode2))
print(readR1.find(barcode2))
print(readR2.rfind(barcode1))
locate1 = readR1.find(barcode1)
locate2 = readR2.rfind(barcode2)
locate3 = readR1.find('KKKKK')
if locate1 < 15 and locate1 != -1:  # locate1 != -1 should be the first check point
    readR1 = readR1[locate1 + 6:]
print(readR1)
if locate2 > 36 and locate2 != -1:  # locate2 != -1 should be the first check point
    readR2 = readR2[:locate2]
print(readR2)
print(readR1[locate3:])

string1 = 'aaaaa\n'
string2 = '...`11!@\r\n'
print(string1.strip(), string1.rstrip(), string1)
print(string2.strip(), string2.rstrip(), string2)
print('aaa')

list1 = ['a', 'c', 'd', 88]
a = 0
list2 = 'aaa'
for x, y in enumerate(list1):
    a += x
    print(x, y)
print(a)
for x, y in enumerate(list2):
    a += x
    print(x, y)
print(a)

g = 'aaa'
a = 0.12345
print('{}\t{:.2f}'.format(g, a))
print('\0', 'aa')
print('\n', 'aa')

print(','.join(map(str, list1)))
import pandas as pd
import numpy as np
import os, sys
df1 = pd.read_table('D://Project/drs/GDC/test.txt', index_col=0)
print(df1.shape)
print(df1.index[0], df1.columns[0])
threshold = 0.2
ii = np.where((df1.values >= threshold) | (df1.values <= -threshold))
print(type(ii))
print(ii)
print(len(ii))
print(len(ii[1]))
for r, c in zip(ii[0], ii[1]):
    print(df1.index[r], df1.columns[c], df1.iloc[r,c])


def rreplace(s, old, new, occurrence):
    li = s.rsplit(old, occurrence)
    return new.join(li)


aaa = 'Corr_pvalue_abcgdgefg'
print(aaa.lstrip('Corr_pvalue_'))   # strip use ''one chr" to trim.
print(rreplace(aaa, 'g', '', 1))

path = '/home/siza/siza/imm/fq_TCR/SRR3233947_1_trim.fastq'
fpath, fname = os.path.split(path)
fbase, fext = os.path.splitext(fname)
srr, others = fname.split('_', 1)
print(srr)
print(path)
print(fpath)
print(fname)
print(fbase)
print(fext)
print(os.path.join(fpath, fbase + '_g' +fext))

string1 = 'TCGA-LUSC_GSE70138_Level4_ZSVCINF_A375.311_m_scipy_pearson_top250_snv_indel.extracted.snp.code.matrix.LUSC_2'
file = string1[string1.rfind('snv_indel'):] + '.csv'
print(file)

df1 = pd.read_table('D:/Project/drs/result_tmp/cmapID2CL.txt', index_col=0)
print(df1.shape)
print(df1)
print(df1.iloc[:2, :])
df1.T.to_csv('D:/Project/drs/result_tmp/cmap_id2clDF', sep='\t')


def openDF(in_path, direct = 'f'):
    fpath, fname = os.path.split(in_path)
    fbase, fext = os.path.splitext(fname)
    df = pd.read_csv(in_path, index_col=0) if fext == '.csv' else pd.read_table(in_path, index_col=0)
    if direct == 't':
        df = df.transpose()
    return df, fbase


pid = 'MCF7'

dfi, dfi_base = openDF('D:/Project/drs/result_tmp/cmap_50.csv')
dfa = pd.read_table('D:/Project/drs/result_tmp/cmap_id2clDF')
dfa2 = pd.read_table('D:/Project/drs/result_tmp/cmapID2CL.txt')
print(dfa)
ii = np.where(dfa.values == pid)[1]
ii2 = np.where(dfa2.values == pid)[0]
print(ii)
print(ii2)
df2 = dfi.iloc[:, ii]
print(df2)


f_1 = 1.3e-6
f_2 = 0.000000001
print(f_1)
print(type(f_1))
print(f_1 - f_2)
print(f_2 + f_1)    # not precision value, 1.3010000000000001e-06

a = False
b = True

if a or b:
    print('1')
else:
    print('2')

import scipy.stats as st

arr1 = np.array([0.3, 0.5, 0.1, 0.12])
arr2 = np.array([0.12, 0.3, 0.1, 0.5, 0.9, 0.92, 0.88, 0.73, 0.48, 0.55])
arr3 = np.array([0.1, 0.0])
t_1 = st.ttest_ind(arr1, arr2)
t_2 = st.ttest_ind(arr2, arr3)
print(t_1)
print(t_2)
print('\t'.join(arr1.astype(str)))
df1 = pd.read_csv('D:/Project/drs/result_tmp/LUSC_A549_250_3_n50.csv', index_col=0)
print('\t'.join(df1.index))

word1 = 'Summary_snv_snv_indel.extracted.snp.LUSC.genelevel.mut.status_2_g_3'
a = word1[word1.find('.snp.') + 5:]
b = a[:a.find('.')]
print(b)

word2 = 'ttest_snv_indel.extracted.snp.LUSC.genelevel.mut.status_2_g_LUSC_A549_m_top100_3.csv'
title = word2.replace('snv_indel.extracted.', '')
title = rreplace(title, '.mut.status', '', 1)
print(title)
gdc = title[title.find('snp.') + 4:]
gdc = gdc[:gdc.find('.')]

print(gdc)

dictA = {'a': 1, 'b': 2, 'c': 3}
dictB = {'d': 4, 'e': 5, 'f': 6}
dictX = {'A': dictA, 'B': dictB}
print(len(dictX), len(dictB), len(dictA))
print(dictX)

string1 = '>aaaa\n'
id = string1.rstrip().split(' ')[0].lstrip('>')
print(id)

string2 = r'C:\read'
print(string2)

list1 = [1, 2, 3, 5, 8]
list2 = [1, 2, 3, 9]
print(set(list1).symmetric_difference(list2))

def cmap_dict(path):
    lineCount = 0
    id2anno = {}
    with open(path) as inputFile:
        for line in inputFile:
            lineCount += 1
            if lineCount == 1:
                continue
            lf = line.rstrip().split('\t')
            drug_annotation = '\t'.join([lf[6], lf[2], lf[13], lf[14], lf[11], lf[4], lf[5]])
            id2anno[lf[0]] = drug_annotation
    return id2anno

id2drug = cmap_dict('D:/Project/drs/gdsc/cmap/cmap_instances_02.txt')
print(len(id2drug))
print(id2drug.get('364'))

list1 = ['aaaTCccc', 'cvdfafa', 'a ferqTCweq']
list2 = [x for x in list1 if 'TC' in x]
print(list2)

string1 = 'aaa\r\n'
string2 = 'bbb\n'
string3 = 'ccc\t'
string4 = 'ddd\r'
print(string3 + string1)
print(string2)
print(string3.rstrip() + string1.rstrip())
print(string2.rstrip())
print(string3.rstrip())
print(string4, string3)

a1 = '14'
a2 = '3'
print('{:.3f}'.format(int(a1) / int(a2)))
l1 = ['aa', '14', 'gg']
l2 = ['gg', 'cc', 'aa']
print('\t'.join(l1 + l2))

q1 = 0.4
q2 = 0.9
q3 = 1.33
print(list(map(str, [q1, q2, q3])))

set1 = set(l1 + l2)
print(set1)
for a in l2:
    if a in set1:
        print(a)

print(int('5') / 3)

dict1 = {'c': 3, 'd': 5, 'a': 0, 'b': 13, 'e': 16, 'h': 'rr', 'g': 1, 'f': '5'}
for k in dict1:
    print('{},{}'.format(k, dict1.get(k)))
for k, v in dict1.items():
    print('{},{}'.format(k, v))
from collections import OrderedDict
dict2 = sorted(dict1)
dict3 = sorted(dict1.items())
dict4 = dict(sorted(dict1.items()))
dict5 = OrderedDict(sorted(dict1.items()))
print(dict1, dict2, dict3, dict4, dict5, sep='\n')

check = '01'
print(int(check))

with open('D:/Project/circ_miRNA/barcode_files/TCGA-ACC.csv') as in_f:
    print(next(in_f))
    #print(in_f.readline().rstrip())
    l_c = 0
    for line in in_f:
        l_c += 1
        if l_c <= 5:
            print(line.rstrip())

list3 = [1, 2, 3, 'a', 5]
list3[3] = 'b'
print(list3)

dict1 = {'g1':{'1': 1, '2': 2, '3':3}, 'g2':{'1':2, '3':'null', '4': 5}, 'g3':{'0': 0, '5': 3}}
print(dict1.keys())
frame1 = pd.DataFrame(dict1)
print(frame1)
frame1 = frame1.replace('null', np.nan)
print(frame1.iloc[:,1])
print(frame1.isnull().sum(1))
print(frame1.T.mean())
frame3 = frame1.T.fillna(frame1.T.mean()).T
print(frame3)


frame1 = frame1.fillna(0)
print(frame1)
print(frame1.sum())
print(frame1.mean(1))
frame2 = (frame1 / frame1.sum()) * 10
print(frame2)
print(frame2.sum())


def open_df(in_path, direct='n'):
    fpath, fname = os.path.split(in_path)
    fbase, fext = os.path.splitext(fname)
    df = pd.read_csv(in_path, index_col=0) if fext == '.csv' else pd.read_table(in_path, index_col=0)
    print('{}: Transpose: {}'.format(fbase, direct))
    if direct == 't':
        df = df.transpose()
    return df, fbase


s1 = pd.Series({'1': 1, '2': 2, '3':3})
s2 = pd.Series({'1':2, '3':'null', '4': 5})
s3 = pd.Series({'0': 0, '5': 3})
df1 = pd.DataFrame()
df1['g1'] = s1
df1['g2'] = s2
print(df1.shape)

str1 = 'AF-141 bST0049 Gis11\n'
print(str1.rstrip().lower().split(' '))

np.random.seed(0)
import matplotlib.pyplot as plt
import csv
df = pd.DataFrame(2.5 * np.random.randn(20,5) + 3)
print(df)
print(df.mean(0), df.std(0))
print(df.iloc[:, 0])

df.iloc[:, 0].plot.density()
plt.show()
str1 = 'aaa'
print(str1)
with open('D:/Project/circ_miRNA/tcga_test.csv') as in_f:
    aaaa, bbbb, cccc, *dddd = csv.reader(in_f)
    print(type(aaaa))
    print(aaaa)
    print(bbbb)
    print(cccc)
    for mir in dddd:
        print(mir)
print(aaaa[0])

import re
target_cancer = 'Prostate'
match_map = {"Prostate": [r"TCGA-PRAD", r"PC\d{1,2}\Z"], "Colon": [r"TCGA-COAD", r"\d{1}S\d{1,2}"], "Healthy": [r"N\d{1,2}\Z"]}
matchs = match_map[target_cancer]
print(matchs)

a = 'atcgnaacg'

if any(c in a for c in ('N', 'n')):
    print('yes')
else:
    print('No')

print(a.replace('aa*', ''))