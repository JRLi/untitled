with open('GSE70138_Level4_ZSVCINF_n115209x22268_20151231_gid.txt') as in_put, open(
        'GSE70138_Level4_ZSVCINF_n115209x22268_20151231_all.txt', 'w') as out_put:
    for line in in_put:
        lf = line.split('\t', maxsplit=1)
        if '///' not in lf[0]:
            out_put.write(line)
        else:
            g_list = lf[0].split('///')
            for gid in g_list:
                out_put.write('\t'.join([gid, lf[1]]))