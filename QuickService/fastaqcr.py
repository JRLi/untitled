#!/usr/bin/env python3.5
"""
Created by JRLi on 2016/12/06 for python learning
"""
import sys
import argparse
import os

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


class FaParser():
    def __init__(self, fa_argv):
        self.fa_argv = fa_argv
        seqFile = self.fa_argv.input

    def test(self):
        print(self.fa_argv)


    def process_seq(self):
        pass


class FqParser():
    def __init__(self, fq_argv):
        self.fq_argv = fq_argv


    def test(self):
        print(self.fq_argv)


def fa_parser(file, base, rev):
    with open('./' + file, 'r') as inputFile:
        id2seq_dict = {}
        id, seq = '', []
        for line in inputFile:
            if line.startswith('>'):
                if id != '':
                    fseq = inverse_complement(''.join(seq), True) if rev else ''.join(seq)
                    fseq = fseq if base == 0 else insert_end(fseq, base)
                    id2seq_dict[id] = fseq
                id = line.replace('\"', '').replace('\n', '')[1:]
                seq = []
            else:
                seq = seq.append(line.replace('\n', ''))
        fseq = inverse_complement(''.join(seq), True) if rev else ''.join(seq)
        id2seq_dict[id] = fseq



def inverse_complement(line, seq = False):
    complement_dict = {'a':'t','t':'a','c':'g','g':'c','A':'T','T':'A','C':'G','G':'C'}
    r_line = line[::-1]
    if seq:
        rc_line = "".join(complement_dict.get(bp, bp) for bp in r_line)
        return rc_line
    else:
        return r_line


def insert_end(line, index):
    return "".join([line[j:j+index]+"\n" for j in range(0, len(line), index)])


def main(argv=None):
    try:
        if argv is None:
            argv = args_parse()
            print('Implement reverse complement.' if argv.reverse else 'No reverse complement.')
            root = '.reverse_complement.' if argv.reverse else ''
            if argv.command is None:
                raise Usage('No command!!')
            elif argv.command == 'fa':
                print('Implement fasta mode.')
                for fileName in argv.input:
                    fbase, fext = os.path.splitext(fileName)
                    with open('./' + fbase + root + fext, 'w') as outFa:
                        id2Seq_dict = fa_parser(fileName, argv.base, argv.reverse)


            else:
                print('Implement fastq mode.')
                fastqc = FqParser(argv)
                fastqc.test()
    except Usage as err:
        print(sys.stderr, err.msg)
        print(sys.stderr, "for help use --help")
        return 2

if __name__ == "__main__":
    sys.exit(main())
