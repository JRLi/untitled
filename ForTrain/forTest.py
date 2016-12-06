# -*- coding: UTF-8 -*-

# import time
# import numpy as np
# from pandas import Series, DataFrame
# from collections import defaultdict
# import pandas as pd
import re
import os
# stringA = 'abcdefg'
# print(stringA[0:-1])
# dictTest = {"a": 1, "b": 2, "c": 3, "e": "aaa"}
# valueTest1 = dictTest.get("d", "")
# valueTest2 = dictTest.get("a", "")
# valueTest3 = dictTest.get("e", "")
# print(valueTest2)
# print(valueTest1)
# print(valueTest3)
# aa = "gi|329755407|gb|HQ832783.1|"
# gg = "gi|329755407|gb|HQ832783.1|HQ832783"
# print(aa.rfind("|"))
# print(aa.rfind("|", 0, -1))
# s1 = aa.rfind("|", 0, -1)
# s2 = aa.rfind("|")
# s3 = gg.rfind("|", 0, s2)
# world1 = gg[gg.rfind("|", 0, gg.rfind("|")) + 1:gg.rfind("|")]
# print(aa[s1 + 1:s2])
# print(aa.rfind("|"), len(aa))
# print(world1)
#
# print()
# bb = ">abcdefgH.2"
# s3 = bb[1:len(bb) - 1]
# print(s3)
# dict2 = {"aa.1": 1, "abc.2": "ss", "cc.2": 3}
# string2 = "abc."
# string5 = "gg"
# string3 = dict2.get(string2, "None")
# string4 = next(v for k, v in dict2.items() if k.startswith(string2))
# string6 = [v for k, v in dict2.items() if k.startswith(string2)]
# if len(string6) == 0:
#     string7 = "-"
# else:
#     string7 = string6[0]
# print(string3)
# print(string2[:-1], string2[1:-1])
# print(string4, string6, len(string6), string7)
#
# testWord = "Septin 5"
# testWord.replace("Septin ", "Spet")
# testWord2 = testWord.replace("Septin ", "Spet")
# print(testWord, testWord2)
# testWord = testWord.replace("Septin ", "Spet")
# print(testWord)
#
# intTest1 = 20
# intTest2 = 4
# print(intTest1 / intTest2, intTest1 // intTest2, intTest1 % intTest2)   # attention, '/' will output float
#
# testWord = "'a'"
# testWord2 = "a"
# testWord4 = "\"a\""     # if the same quotation (all single or all double), need to use '\' marks
# testWord3 = """b
# a"""                # even print it with \n symbol
# print(testWord, testWord4, testWord2, testWord3)
#
# testWord = "a,b,c,d,e,f,g"
# testField = testWord.split(",")     # important usage
# testWord2 = ",".join(testField)     # important usage, use 'a common string' to combine 'list' to a 'long string'
# print(testField)
# print(testWord2)
# print(testWord.replace(",", ""))
# print(testWord.replace(",", "", 2))     # important, can use this method to replace partial symbol!!!
# del testWord, testWord2     # testWord and testWord2 will be not defined
#
# print("dict a")
# dictA = {"a": [11, 3, 5], "c": [2], "e": [3, 4]}
# listA = dictA.get("b", [])
# listA.append(6)
# dictA["b"] = listA
# print(dictA)
# listA = dictA.get("b", [])
# listA.append(8)
# dictA["b"] = listA
# print(dictA)
# outList = []
# for line in outList:
#     print(line)
#
# testList = ['a', 'b', 'c', 'd']
# for i in range(len(testList)):      # range can omit the start
#     print(testList[i])
#
# tStart = time.time()    # 計時開始
# print(tStart)
# time.sleep(1)
# print("abc")
# for x in range(1000):
#     x += 1
# print(x)
#
# tEnd = time.time()  # 計時結束
#
# print("It cost %f sec" % (tEnd - tStart))   # 會自動做近位
# print(tEnd - tStart)    # 原型長這樣
#
#
# data1 = [1, 7.5, 8, 0, 6, 1]
# print("index of 1", data1.index(1, 0, 1))
# arr1 = np.array(data1)
# print(arr1)
# obj = Series([4, -5, 7, 3])
# print(obj)
# print()
# orgID = 'hsa'
# hrefID = '<a href="/dbget-bin/www_bget?hsa:130589">130589</a>  a href="/dbget-bin/www_bget?hsa:4422">4422</a>'
# print(hrefID)
# pat = '\/dbget-bin\/www_bget\?' + orgID + ':\d+'
# match = re.findall(pat, hrefID)
# if match:
#     for stringG in match:
#         print(stringG[stringG.rfind(":") + 1:])
# else:
#     print("no match")
#
# print()
#
# with open("D:/Project/orchid/KEGGpyTest/pythonReadlinesTest") as inputFile:
#     line = inputFile.read().replace("\n", "")
#     print(inputFile, "type:", type(inputFile))          # read()後便不能readlines(), 反之亦然
#     print(inputFile.readlines(), "type:", type(inputFile.readlines()))
#
#     print(line, "type:", type(line))
#
# stringT1 = 'test'
# stringT2 = 'testaaa'
# print(stringT1 == stringT2)
# print(stringT1 != stringT2)
# print(stringT1 != stringT2[:-3])
# print(stringT1 == stringT2[:-3])
#
# xxx = 10
#
#
# def func():
#     xxx = 15
#     print(xxx)
#     return xxx
#
# yyy = func()
# print(yyy)
# xdict = {}
# print(xdict.keys())
# a = xdict.get('p63', 'x')
# print(a)
# print(xdict.keys())
# ydict = defaultdict(int)
# aa = ydict['p73']
# print(ydict.keys(), aa)
# geneid2godict = defaultdict(lambda: set())
# print(geneid2godict.keys())
# print(type(geneid2godict['p53']))
# print(geneid2godict.keys())
# a = geneid2godict['p53']
# a.add("14")
# a.add("20")
# print(geneid2godict['p53'])
# print(geneid2godict.keys())

mdict = {}
adict = {}
zdict = {'a': 'fff'}
pdict = {'b': 'zzz'}
x = 15
y = 15


def test():
    x = 20
    global y, pdict
    y += 5
    sdict = dict()  # dictionary declared in the function, but can be return!!
    sdict['sss'] = 'aaa'
    mdict['a'] = 'ccc'
    mdict['b'] = 'eee'
    adict['a'] = 'ggg'
    zdict = {'b': 'vvv'}
    pdict = {'g': 'rrr'}
    print('local function x:', x)
    print('local dict zdict:', zdict)
    return mdict, sdict
aaa = test()

print(mdict)
print(adict)
print('not use global, assignment in no use out of function:', zdict)
print('use global, can be assignment in function:', pdict)
print(x)
print(y)
print('because return 2 dictionary:', aaa)

aaaa = ''
bbbb = None
cccc = 'xxx'
print(aaaa, bbbb)
if aaaa:
    print("yes aaaa")
else:
    print("no aaaa")

print('yes bbbb' if bbbb else 'no bbbb')
print(bool(aaaa), bool(bbbb), bool(cccc))

aaa = range(10, 20, 2)
for a in aaa:
    print(a)

stringTest = 'file_record_transcript.pdfxxx1e123 file_07241999.pdf 1231'
pattern = '(file_\w+)(?:(?:t\.pdf)|(?:9\.pdf))'
pattern2 = 'pdf$'
match = re.search(pattern, stringTest)
match2 = re.findall(pattern, stringTest)
if match is None:
    print('none')
elif not match:
    print('no result')
else:
    print(match.group())

print(match2)
for line in match2:
    print(line)

tuple1 = 1, 'aa', 33, 'bbb'
b1, b2, b3 = 1, 2, '333'
a1 = 1; a2 = 2; a3 = '333'
# tuple1[2] = 3   # TypeError: 'tuple' object does not support item assignment
print(tuple1[2])
b1 = 5
print(b1)   # assignment ok
a1 = tuple1[2]
print(a1, type(a1))

line = '.abcde'
fx = line.find('.')
fy = line.find('cde')
fz = line.find('g')
print(fx, fy, fz)

line = 'aaaa\tbbbb\tcccc\n'
item1, item2 = line.strip().split('\t', 1)      # important
print(item1, item2, sep="\n")

filepath = 'http://rest.kegg.jp/list/pathway/test1.fastq.gz'
filepath2 = './gene_pathway_out/test2_aaa_a1.3302.fq.gz'
filepath3 = 'aea_5_r2.fa.gz'
fpath, fname = os.path.split(filepath)
print(fpath, fname, sep="\n")
fbase, fext = os.path.splitext(fname)
print(fbase, fext, sep='\n')

fpath, fname = os.path.split(filepath2)
print(fpath, fname, sep="\n")
fbase, fext = os.path.splitext(fname)
print(fbase, fext, sep='\n')

fpath, fname = os.path.split(filepath3)
print(fpath, fname, sep="\n")
fbase, fext = os.path.splitext(fname)
print(fbase, fext, sep='\n')
