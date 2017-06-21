#!/usr/bin/env python
with open('CCLE_hybrid_capture1650_hg19_allVariants_2012.05.07.maf') as input_file, open('allv_cDNAchange','w') as outFile:
    for line in input_file:
        lf = line.replace('\r\n', '').split('\t')
        if lf[35] != '':
            outFile.write(line)