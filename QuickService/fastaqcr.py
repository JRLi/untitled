#!/usr/bin/env python3.5
"""
Created by JRLi on 2016/12/06 for python learning
"""
import sys
import getopt
import argparse
import gzip
import os

use_message = '''
    To calculate NGS reads length and the number of reads.
    To perform the reverse complement of fasta/fastq file.
    Uniform fasta file line length: one line or 60mer/per line.
'''


def args_parse():
    parser = argparse.ArgumentParser(description=use_message, epilog="  Attention: For test!!",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-r", "--reverse", action="store_true", help="output reverse complement")
    subparsers = parser.add_subparsers(help='choice file format')
    parser_a = subparsers.add_parser('fa', help='input is fasta format')
    parser_a.add_argument('-b', '--base', type=int, default=60, help='base pairs per line; if set 0, all in one line')
    parser_a.add_argument('input', nargs='+', help="input fasta file")
    parser_b = subparsers.add_parser('fq', help='input is fastq format')
    parser_b.add_argument('-c', choices='xyz', help='baz help')

    args = parser.parse_args()
    return args


class Usage(Exception):
    def __init__(self, msg):
        print()
        self.msg = msg


class Params:
    def __init__(self, type):
        self.type = type

    def parse_options(self, argv):
        try:
            opts, args = getopt.getopt(argv[1:], "ho:", ["help", "output="])

        except getopt.error as msg:
            raise Usage(msg)
        return opts, args


def main(argv=None):
    try:
        if argv is None:
            argv = args_parse()
            print('Implement reverse complement' if argv.reverse else 'Normal mode')
            root = '.reverse_complement.' if argv.reverse else ''
        print(argv)

    except Usage as err:
        print(sys.stderr, err.msg)
        print(sys.stderr, "for help use --help")
        return 2

if __name__ == "__main__":
    sys.exit(main())
