#!/usr/bin/env python3.5
"""
Created by JRLi on 2016/12/06 for python learning
"""
import sys
import argparse
import os
from collections import OrderedDict

use_message = '''
    To calculate NGS reads length and the number of reads.
    To perform the reverse complement of fasta/fastq file.
    Uniform fasta file line length: one line or 60mer/per line.
'''


def args_parse():
    parser = argparse.ArgumentParser(prog='fastaqcr', description=use_message, epilog="  Attention: For test!!",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-r", "--reverse", action="store_true", help="output reverse complement")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s beta2')
    subparsers = parser.add_subparsers(help='choice file format', dest='command')

    parser_a = subparsers.add_parser('fa', help='input is fasta format')
    # parser_a.add_argument('-b', '--base', type=int, default=60, help='base pairs per line; if set 0, all in one line')
    parser_a.add_argument('-b', '--base', type=int, help='base pairs per line; if set 0, all in one line')
    parser_a.add_argument('-i', '--id2seq', action="store_true", help="output ID2Sequence file, separated by \\t")
    parser_a.add_argument('input', nargs='+', help="input fasta file")
    parser_a.add_argument('-s', '--split', type=int, default=1, help='split fasta file two [-s] parts')
    parser_b = subparsers.add_parser('fq', help='input is fastq format')
    parser_b.add_argument('-p', dest='paired', nargs=2,help="Comma-separated list of files containing the #1 mates, "
                                                            "and space-separated list of paired #2 mates, "
                                                            "E.g, x_1.fq,y_1.fq,z_1.fa x_2.fq,y_2.fq,z_2.fa",
                          metavar=('Left', 'Right'))
    # parser_b.add_argument('-1', dest='left', nargs='?', help="Comma-separated list of files containing the #1 mates, "
    #                                                      "no space inside. E.g, x_1.fq,y_1.fq,z_1.fa")
    # parser_b.add_argument('-2', dest='right', nargs='?', help="Comma-separated list of files containing the #2 mates, "
    #                                                        "no space inside. .g, x_2.fq,y_2.fq,z_2.fa")
    parser_b.add_argument('-s', dest='single', nargs='?', help='Comma-separated list of single-end files')
    args = parser.parse_args()
    return args


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


def fa_parser(file, base, rev):
    with open(file, 'r') as inputFile:
        id2seq_dict = OrderedDict()
        id, seq, total_base, total_gc, len_list, check_len, n50, n90 = '', [], 0, 0, [], 0, [], []
        for line in inputFile:
            if line.startswith('>'):
                if id != '':
                    fseq = reverse_complement(''.join(seq), True) if rev else ''.join(seq)
                    total_base += len(fseq)
                    len_list.append(len(fseq))
                    total_gc += gc_base(fseq)
                    bfseq = fseq + '\n' if base in (0, None) else insert_end(fseq, base)
                    id2seq_dict[id] = bfseq
                id = line.replace('\"', '').replace('\r\n','').replace('\n', '')[1:]
                seq = []
            else:
                seq.append(line.replace('\r\n','').replace('\n', ''))
        fseq = reverse_complement(''.join(seq), True) if rev else ''.join(seq)
        total_base += len(fseq)
        len_list.append(len(fseq))
        total_gc += gc_base(fseq)
        bfseq = fseq + '\n' if base in (0, None) else insert_end(fseq, base)
        id2seq_dict[id] = bfseq
        len_list.sort(reverse=True)
        max_contig, min_contig = len_list[0], len_list[-1]
        gc_content = total_gc/total_base
        total_contig = len(id2seq_dict)
        for j in len_list:
            check_len += j
            if check_len >= total_base/2:
                n50.append(j)
            if check_len >= total_base*9/10:
                n90.append(j)
        return id2seq_dict, total_contig, max_contig, min_contig, n50[0], n90[0], total_base, gc_content


def dict2file(dictIn, arg_id2seq, base_number, attlist):
    if arg_id2seq:
        with open(attlist[0] + attlist[1] + attlist[3] + '.id2seq', 'w') as inputFile:
            for k, v in dictIn.items():
                inputFile.write(k + '\t' + v.replace('\r\n','').replace('\n', '') + '\n')
    if base_number is not None:
        with open(attlist[0] + attlist[1] + '.' + str(base_number) + attlist[3] + attlist[2], 'w') as outFa:
            for id in dictIn.keys():
                outFa.write('>' + id + '\n' + dictIn[id])


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


def gc_base(line):
    return len(''.join(bp for bp in line if bp in ['c', 'g', 'C', 'G']))


def fq_processor(path, rev, file_in):
    ReadCount, lineTemp, total_base, total_gc = 0, 0, 0, 0
    with open(path + file_in, 'r') as readFileL:
        if rev:
            fbase, fext = os.path.splitext(file_in)
            outL = open(path + fbase + '.reverse_complement' + fext, 'w')
        for line in readFileL:
            lineTemp += 1
            if lineTemp in {1, 3}:
                if rev:
                    outL.write(line)
            if lineTemp == 2:
                seq_temp = line.replace('\n', '').replace('\"', '').strip()
                total_base += len(seq_temp)
                total_gc += gc_base(seq_temp)
                ReadCount += 1
                if rev:
                    outL.write(reverse_complement(seq_temp, True) + '\n')
            if lineTemp ==4:
                if rev:
                    outL.write(reverse_complement(line.replace('\n', '').strip()) + '\n')
                lineTemp = 0
        average_length = total_base/ReadCount
        gc_content = total_gc/total_base
        if rev:
            outL.close()
    return ReadCount, average_length, total_base, gc_content


def main(argv=None):
    try:
        if argv is None:
            argv = args_parse()
            print('Implement reverse complement.' if argv.reverse else 'No reverse complement.')
            root = '.reverse_complement' if argv.reverse else ''
            if argv.command is None:
                raise Usage('No command!!\nMust input fa or fq!')
            elif argv.command == 'fa':
                print('Implement fasta mode.')
                print(argv)
                with open('./SummaryOfFasta.txt', 'w') as rs:
                    for fileName in argv.input:
                        fpath, fname = os.path.split(fileName)
                        fbase, fext = os.path.splitext(fname)
                        id2Seq_dict, tc, mxc, mic, n50, n90, tb, gcc = fa_parser(fileName, argv.base, argv.reverse)
                        if argv.id2seq is not None or argv.base is not None:
                            dict2file(id2Seq_dict, argv.id2seq, argv.base, [fpath, fbase, fext, root])
                        rs.write('File: {}\nTotal_contigs: {}\nMax_contig: {}\nMin_contig: {}\nn50: {}\nn90: {}\n'
                                 'Total_bases: {}\nGC_content: {}\n\n'.format(fname, tc, mxc, mic, n50, n90, tb, gcc))
            else:
                print('Implement fastq mode.')
                print(argv)
                if argv.paired is None and argv.single is None:
                    raise Usage('Error! No input file, please use:./fastaqcr [-r] fq [-p PAIRED PAIRED] [-s [SINGLE]]')
                with open('./SummaryOfReads.txt', 'w') as results:
                    if argv.single is not None:
                        print('Run single end reads.\nInput line: {}\n'.format(argv.single))
                        for sfile in argv.single.split(','):
                            fpath, fname = os.path.split(sfile)
                            results.write(fname + '\t' + 'Single end\n')
                            rc, al, tb, gc = fq_processor(fpath, argv.reverse, fname)
                            results.write('Reads_count:{}\nAverage_read_length:{}\nTotal_bases:{}\nGC_content:{}\n\n'.
                                          format(rc, al, tb, gc))
                    if argv.paired is not None:
                        print('Run paired end reads.\nInput line: {}\n'.format(argv.paired))
                        r1_list = argv.paired[0].split(',')
                        r2_list = argv.paired[1].split(',')
                        if len(r1_list) != len(r2_list):
                            print('Warning: the number of r1 files is not equal to r2 files!!!\n'
                                  'Warning: the definition of paired-end is according to shorter list!!!')
                        for left, right in zip(r1_list, r2_list):
                            fpathL, fnameL = os.path.split(left)
                            fpathR, fnameR = os.path.split(right)
                            print('Processing paired:\t{}')
                            results.write(fnameL + '\t' + 'Paired end R1\n')
                            rc_r1, al_r1, tb_r1, gc_r1 = fq_processor(fpathL, argv.reverse, fnameL)
                            rc_r2, al_r2, tb_r2, gc_r2 = fq_processor(fpathR, argv.reverse, fnameR)
                            results.write('Reads_count:{}\nAverage_read_length:{}\nTotal_bases:{}\nGC_content:{}\n'.
                                          format(rc_r1, al_r1, tb_r1, gc_r1))
                            results.write(fnameR + '\t' + 'Paired end R2\n')
                            results.write('Reads_count:{}\nAverage_read_length:{}\nTotal_bases:{}\nGC_content:{}\n\n'.
                                          format(rc_r2, al_r2, tb_r2, gc_r2))

    except Usage as err:
        print(sys.stderr, err.msg)
        print(sys.stderr, "Terminated, for help use -h or --help")
        return 2

if __name__ == "__main__":
    sys.exit(main())
