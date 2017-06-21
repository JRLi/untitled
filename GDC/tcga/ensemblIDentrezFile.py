#!/usr/bin/env python3.6
from collections import defaultdict

hgncPath = 'hgnc_id_170512'
emartPath = 'ensembl_eid2gid_mart_export.txt'
gdcData = 'rna_seq_sampleID_FPKM.csv'


def main():
    ensembl2entrez_dict = defaultdict(set)
    with open(hgncPath) as hgnc_in, open(emartPath) as emart_in:
        for line in emart_in:
            if line.startswith('Gene'):
                continue
            lf = line.replace('\r\n', '').replace('\n', '').replace('\"', '').split('\t')
            if lf[2] != '':
                ensembl2entrez_dict[lf[0]].add(lf[2])
        for line in hgnc_in:
            if line.startswith('hgnc_id'):
                continue
            lf = line.replace('\r\n', '').replace('\n', '').replace('\"', '').split('\t')
            if lf[4] != '' and lf[5] != '':
                ensembl2entrez_dict[lf[5]].add(lf[4])
        print(len(ensembl2entrez_dict))
        ambiguous = 0
        for k, v in ensembl2entrez_dict.items():
            if len(v) > 1:
                ambiguous += 1
                print(k, v)
        print(ambiguous)
    line_count, noID_count, entrez_count = 0, 0, 0
    with open(gdcData) as gdc_in, open('rna_seq_gID_sampleID_FPKM.csv', 'w') as gdc_out, open('noEid', 'w') as no_out:
        for line in gdc_in:
            line_count += 1
            if line_count == 1:
                gdc_out.write(line)
                continue
            ensemblPre, exp_value = line.split(',', maxsplit=1)
            ensembl_id = ensemblPre[:ensemblPre.rfind('.')]
            entrez_set = ensembl2entrez_dict.get(ensembl_id, '')
            if entrez_set == '':
                noID_count += 1
                no_out.write(ensembl_id + '\n')
            else:
                for entrez_id in entrez_set:
                    entrez_count += 1
                    gdc_out.write(','.join([entrez_id, exp_value]))
    print('line, noID, ent_id:', line_count, noID_count, entrez_count)

if __name__ == '__main__':
    main()