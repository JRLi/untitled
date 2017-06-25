#!/usr/bin/env python3.6
import argparse
import os
import sys


def args_parse():
    parser = argparse.ArgumentParser(description='for barcode trimming')
    parser.add_argument('-f', '--srr', type=str, default='srr_list', help='srr 2 title')
    parser.add_argument('-b', '--barcode', type=str, default='barcode_list.txt', help='barcode 2 seq')
    parser.add_argument('-s', dest='single', nargs='?', help='A single-end files')
    parser.add_argument('-p', dest='paired', nargs=2,help="file_1.fq file_2.fq", metavar=('R1', 'R2'))
    parser.add_argument('-l', '--thresholdL', type=int, default=12, help="barcode before 5' nt threshold")
    parser.add_argument('-r', '--thresholdR', type=int, default=36, help="barcode_rc behind 3' nt threshold")
    args = parser.parse_args()
    return args


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


def srr2barcode_dict(barcode_annotation, SRR_annotation):
    srr2seq_dict = {}
    with open(barcode_annotation) as barcodeFile, open(SRR_annotation) as srrFile:
        bc2seq_dict = dict(line.strip().split('\t', 1) for line in barcodeFile)
        for line in srrFile:
            lf = line.strip().split('\t')
            bc = lf[1].split('.')[1] if lf[1].startswith('TCR') else ''
            srr2seq_dict[lf[0]] = bc2seq_dict.get(bc, '')
    return srr2seq_dict


def reverse_complement(line, seq=False):
    complement_dict = {'a':'t','t':'a','c':'g','g':'c','A':'T','T':'A','C':'G','G':'C'}
    r_line = line[::-1]
    if seq:
        rc_line = "".join(complement_dict.get(bp, bp) for bp in r_line)
        return rc_line
    else:
        return r_line


def barcode_trimmer(path, barcode, output, section, thresholdL, thresholdR):
    ReadCount, lineTemp, processed, locate = 0, 0, 0, 0
    with open(path) as readFile, open(output, 'w') as trimFile:
        for line in readFile:
            lineTemp += 1
            if lineTemp in [1, 3]:
                trimFile.write(line)
            if lineTemp == 2:
                ReadCount += 1
                locate = line.strip().find(barcode) if section == 0 else line.strip().rfind(reverse_complement(barcode, True))
                if section == 0:
                    if locate <= thresholdL and locate != -1:
                        processed += 1
                        read = line.strip()[locate + 6:]
                        trimFile.write(read + '\n')
                    else:
                        trimFile.write(line)
                else:
                    if locate >= thresholdR and locate != -1:
                        processed += 1
                        read = line.strip()[:locate]
                        trimFile.write(read + '\n')
                    else:
                        trimFile.write(line)
            if lineTemp == 4:
                if section == 0:
                    if locate <= thresholdL and locate != -1:
                        read = line.strip()[locate + 6:]
                        trimFile.write(read + '\n')
                    else:
                        trimFile.write(line)
                else:
                    if locate >= thresholdR and locate != -1:
                        read = line.strip()[:locate]
                        trimFile.write(read + '\n')
                    else:
                        trimFile.write(line)
                lineTemp = 0

    return ReadCount, processed


def path_parser(path, index, srr2seq_dict, five_prime_th, three_prime_th):
    fpath, fname = os.path.split(path)
    srr, others = fname.split('_', 1)
    print(srr)
    outfile = '_'.join([srr, 'tb', others])
    barcode = srr2seq_dict.get(srr, '')
    print(barcode)
    if barcode != '':
        reads, process = barcode_trimmer(path, barcode, outfile, index, five_prime_th, three_prime_th)
        print('{}\t{}\tall:{}\ttb:{}\trate:{:.3f}'.format(outfile, barcode, reads, process, process / reads))
    else:
        print(srr + ' has no barcode.')


def main(argv=None):
    if argv is None:
        argv = args_parse()
        print(argv)
        srr2seq_dict = srr2barcode_dict(argv.barcode, argv.srr)
        print(len(srr2seq_dict))
        if argv.paired is not None:
            for index, file_path in enumerate(argv.paired):
                path_parser(file_path, index, srr2seq_dict, argv.thresholdL, argv.thresholdR)
        if argv.single is not None:
            path_parser(argv.single, 0, srr2seq_dict, argv.thresholdL, argv.thresholdR)

if __name__ == "__main__":
    sys.exit(main())