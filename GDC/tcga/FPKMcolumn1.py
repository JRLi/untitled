#!/usr/bin/env python3.6


def uuid2barcode_col1():
    uuid2caseDict, uuid2sampleDict = {}, {}
    with open('all_37_proj_fpkm.csv') as in_file:
        for line in in_file:
            lf = line.split(',')
            uuid2caseDict[lf[0]] = lf[2]
            sample = '-'.join(lf[2].split('-')[0:4])
            uuid2sampleDict[lf[0]] = sample
    outList, outList2 = [], []
    with open('col1.csv') as in_file, open('mod_col1.csv', 'w') as out_file, open('s_col1.csv', 'w') as out_file2:
        for line in in_file:
            lf = line.replace('\r\n', '').replace('\n', '').split(',')
            for index, uuid in enumerate(lf):
                caseID = uuid2caseDict.get(uuid, '')
                sampleID = uuid2sampleDict.get(uuid, '')
                outList.append(caseID)
                outList2.append(sampleID)
                if caseID == '' or sampleID == '':
                    print(index, uuid)
        out_file.write(','.join(outList) + '\n')
        out_file2.write(','.join(outList2) + '\n')


def case2sample_col1():
    with open('newcol1.csv') as in_file, open('s_col1.csv', 'w') as out_file:
        outList = []
        for line in in_file:
            lf = line.replace('\r\n', '').replace('\n', '').split(',')
            for caseID in lf:
                if caseID == '':
                    outList.append(caseID)
                else:
                    cl = caseID.split('-')
                    sampleID = '-'.join(cl[0:4])
                    outList.append(sampleID)
        print(len(outList))
        out_file.write(','.join(outList) + '\n')


def s2p_from_c2p():
    with open('c2p.csv') as in_file, open('s2p.csv', 'w') as out_file:
        for line in in_file:
            cid, pid = line.replace('\r\n', '').replace('\n', '').split(',')
            sid = '-'.join(cid.split('-')[0:4])
            out_file.write(','.join([sid, pid]) + '\n')


def s2p_final():
    with open('s2p.csv') as s2pFile, open('sample_id') as sidFile, open('s2pDF.txt', 'w') as outFile:
        s2p_dict = {}
        p_list = []
        s_list = []
        for line in s2pFile:
            lf = line.replace('\n', '').split(',')
            s2p_dict[lf[0]] = lf[1]
        for line in sidFile:
            sid = line.replace('\r\n', '').replace('\n', '')
            s_list.append(sid)
            pid = s2p_dict.get(sid, '')
            p_list.append(pid)
            if pid == '':
                print('no pid:', sid)
        outFile.write('\t'.join(s_list) + '\n')
        outFile.write('\t'.join(p_list) + '\n')


def match2linc_gid():
    with open('Lincs_gid') as lincgid_in, open('rna_seq_gID_sampleID_FPKM_m.csv') as gdc_in, \
            open('rna_seq_commonID_sampleID_FPKM_m.csv', 'w') as gdc_out:
        lincs_gid_list = [x.replace('\n', '') for x in lincgid_in]
        line_count, noID_count, withID_count = 0, 0, 0
        for line in gdc_in:
            line_count += 1
            if line_count == 1:
                gdc_out.write(line)
                continue
            gid, exp_value = line.split(',', maxsplit=1)
            if gid in lincs_gid_list:
                withID_count += 1
                gdc_out.write(line)
            else:
                noID_count += 1
        print('line:{}\tcommonID:{}\tnoID:{}\n'.format(line_count, withID_count, noID_count))


def main():
    #case2sample_col1()
    #match2linc_gid()
    #s2p_from_c2p()
    s2p_final()

if __name__ == '__main__':
    main()