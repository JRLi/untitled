#!/usr/bin/env python3.5
"""
Created by JRLi on 2016/12/03 for python learning
"""
import argparse
import sys

use_message = '''
    Just for arg_parse test
'''


def args_parse():
    parser = argparse.ArgumentParser(description=use_message)
    parser.add_argument("-r", "--reversCP", action="store_true", help="output revers complement")
    # parser.add_argument('organism', nargs='+', help="kegg org code")
    parser.add_argument("-v", "--verbose", action="store_true", help="output verbose results")
    subparsers = parser.add_subparsers(help='choice file format')
    parser_a = subparsers.add_parser('fa', help='input is fasta format')
    parser_a.add_argument('organism', nargs='+', help="kegg org code")
    parser_b = subparsers.add_parser('fq', help='input is fastq format')
    parser_b.add_argument('-b', choices='xyz', help='baz help')

    args = parser.parse_args()
    return args


def main(argv=None):
    if argv is None:
        argv = args_parse()
        print(argv)
        # print('Number of organism:', len(argv.organism), '\n', 'Detail:', argv.organism)
        print('Verbose mode' if argv.verbose else 'Default quiet mode')

if __name__ == "__main__":
    sys.exit(main())
