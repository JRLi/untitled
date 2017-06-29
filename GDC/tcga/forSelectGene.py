#!/usr/bin/env python
import os
import sys

delimiter = ','
outpath = '.'
in_path = '../snv_gene_2'


def getGeneList(path):
    with open(path) as file_in:
        gene_list = [x.split(delimiter)[0] for x in file_in if not x.startswith('Gene Symbol')]
    return gene_list


def out_filter(path, g_list):
    fpath, fname = os.path.split(path)
    fbase, fext = os.path.splitext(fname)
    with open(path) as file_in, open(os.path.join(outpath, fbase + '_g' +fext), 'w') as out_file:
        lc = 0
        for line in file_in:
            if line.startswith(delimiter):
                out_file.write(line)
            else:
                gene_check = line.split(delimiter, 1)[0].split(':')[0]
                if gene_check in g_list:
                    out_file.write(line)
                    lc += 1
        print(lc)


def main():
    geneList = getGeneList('Census_allThu_Jun29_2017.csv')
    print(len(geneList))
    file_list = os.listdir(in_path)
    for file in file_list:
        if file.startswith('snv_indel'):
            out_filter(os.path.join(in_path, file), geneList)

if __name__ == '__main__':
    main()