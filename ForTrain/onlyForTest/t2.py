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
seq_ill_uni =   'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'
seq_3 = 'CGCTCTTGGGAACACGT' # SRR3234010_1.fastq have 9711954 at 4st nt
seq_2_pp2 = 'AGAGGACCTGAACAACGTGTTCCCAAGA'
read1='AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'
read2 = 'TGCAGGGCCAGCA'
seqr = reverse_complement(read2, True)
print(seqr)

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
