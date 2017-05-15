#!/usr/bin/env python
from collections import defaultdict
import sys
import os

def main():
    exp_file = sys.argv[1]
    fpath, fname = os.path.split(exp_file)
    fbase, fext = os.path.splitext(fname)

    gene2cell2exp_dict = defaultdict(defaultdict)
    cell_set, gene_set = set(), set()

    with open(exp_file) as inFile, open(fbase + '_df.tsv','w') as outFile:
        for line in inFile:
            c_id, g_symbol, exp = line.replace('\r\n', '').replace('\n', '').split('\t')
            gene2cell2exp_dict[g_symbol][c_id] = exp
            cell_set.add(c_id)
            gene_set.add(g_symbol)
        print('cell number:', len(cell_set))
        print('gene number:', len(gene2cell2exp_dict))

        cell_list, gene_list = list(cell_set), list(gene_set)
        cell_list.sort()
        gene_list.sort()
        outFile.write('\t' + '\t'.join(cell_list) + '\n')

        for geneSymbol in gene_list:
            exp_list = [gene2cell2exp_dict[geneSymbol][c] for c in cell_list]
            outFile.write(geneSymbol + '\t' + '\t'.join(exp_list) + '\n')

if __name__ == '__main__':
    main()