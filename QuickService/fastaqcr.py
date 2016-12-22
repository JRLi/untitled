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
    parser.add_argument('-v', '--version', action='version', version='%(prog)s beta1')
    subparsers = parser.add_subparsers(help='choice file format', dest='command')

    parser_a = subparsers.add_parser('fa', help='input is fasta format')
    parser_a.add_argument('-b', '--base', type=int, default=60, help='base pairs per line; if set 0, all in one line')
    parser_a.add_argument('-i', '--id2seq', action="store_true", help="output ID2Sequence file, separated by \\t")
    parser_a.add_argument('input', nargs='+', help="input fasta file")
    parser_a.add_argument('-s', '--split', type=int, default=1, help='split fasta file two [-s] parts')
    parser_b = subparsers.add_parser('fq', help='input is fastq format')
    parser_b.add_argument('-1', dest='left', nargs='?', help="Comma-separated list of files containing the #1 mates, "
                                                         "no space inside. E.g, x_1.fq,y_1.fq,z_1.fa")
    parser_b.add_argument('-2', dest='right', nargs='?', help="Comma-separated list of files containing the #2 mates, "
                                                           "no space inside. .g, x_2.fq,y_2.fq,z_2.fa")
    args = parser.parse_args()
    return args


class Usage(Exception):
    def __init__(self, msg):
        print()
        self.msg = msg


def fa_parser(file, base, rev):
    with open('./' + file, 'r') as inputFile:
        id2seq_dict = OrderedDict()
        id, seq = '', []

        for line in inputFile:
            if line.startswith('>'):
                if id != '':
                    fseq = reverse_complement(''.join(seq), True) if rev else ''.join(seq)
                    bfseq = fseq + '\n' if base == 0 else insert_end(fseq, base)
                    id2seq_dict[id] = bfseq
                id = line.replace('\"', '').replace('\r\n','').replace('\n', '')[1:]
                print(id)
                seq = []
            else:
                seq.append(line.replace('\r\n','').replace('\n', ''))

        fseq = reverse_complement(''.join(seq), True) if rev else ''.join(seq)
        bfseq = fseq + '\n' if base == 0 else insert_end(fseq, base)
        id2seq_dict[id] = bfseq
        return id2seq_dict


def dict2file(dictIn, pathAndName):
    with open(pathAndName, 'w') as inputFile:
        for k, v in dictIn.items():
            inputFile.write(k + '\t' + v.replace('\r\n','').replace('\n', '') + '\n')


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


def fq_processor(path, rev, file_r1):
    ReadCount, lineTemp, total_base, total_gc = 0, 0, 0, 0
    with open(path + file_r1, 'r') as readFileL:
        if rev:
            global root
            fbase_r1, fext_r1 = os.path.splitext(file_r1)
            outL = open(path + fbase_r1 + root + fext_r1, 'w')
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
                raise Usage('No command!!')

            elif argv.command == 'fa':
                print('Implement fasta mode.')
                print(argv)
                for fileName in argv.input:
                    fbase, fext = os.path.splitext(fileName)
                    with open('./' + fbase + '_' + str(argv.base) + root + fext, 'w') as outFa:
                        id2Seq_dict = fa_parser(fileName, argv.base, argv.reverse)
                        if argv.id2seq:
                            dict2file(id2Seq_dict, './' + fbase + '_id2seq')
                        for id in id2Seq_dict.keys():
                            outFa.write('>' + id + '\n' + id2Seq_dict[id])

            else:
                print('Implement fastq mode.')
                print(argv)
                with open('./SummaryOfReads.txt', 'w') as results:
                    if argv.right is None:
                        print('Run single end reads.')
                        for single in argv.left.split(','):
                            fpath, fname = os.path.split(single)
                            results.write(fname + '\n' + 'Single end\n')
                            rc, al, tb, gc = fq_processor(fpath, argv.reverse, fname)
                            results.write('Reads_count:{}\nAverage_read_length:{}\nTotal_bases:{}\nGC_content:{}\n\n'.
                                          format(rc, al, tb, gc))
                    else:
                        print('Run paired end reads.')
                        for left, right in zip(argv.left.split(','), argv.right.split(',')):
                            fpathL, fnameL = os.path.split(left)
                            fpathR, fnameR = os.path.split(right)
                            results.write(fnameL + '\t' + 'Paired end R1\n')
                            rc_r1, al_r1, tb_r1, gc_r1 = fq_processor(fpathL, argv.reverse, fnameL)
                            rc_r2, al_r2, tb_r2, gc_r2 = fq_processor(fpathR, argv.reverse, fnameR)
                            results.write('Reads_count:{}\nAverage_read_length:{}\nTotal_bases:{}\nGC_content:{}\n'.
                                          format(rc_r1, al_r1, tb_r1, gc_r1))
                            results.write(fnameR + '\t' + 'Paired end R2\n')
                            results.write('Reads_count:{}\nAverage_read_length:{}\nTotal_bases:{}\nGC_content:{}\n'.
                                          format(rc_r2, al_r2, tb_r2, gc_r2))

    except Usage as err:
        print(sys.stderr, err.msg)
        print(sys.stderr, "for help use --help")
        return 2

if __name__ == "__main__":
    sys.exit(main())
