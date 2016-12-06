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


'''


def args_parse():
    parser = argparse.ArgumentParser(description=use_message)
    print("parser after argparse.ArgumentParser:", parser, "type", type(parser))
    parser.add_argument("-r", "--reversCP", action="store_true", help="output revers complement")
    parser.add_argument("-v", "--verbose", action="store_true", help="output verbose results")
    subparsers = parser.add_subparsers(help='choice file format')
    parser_a = subparsers.add_parser('fa', help='input is fasta format')
    parser_a.add_argument('organism', nargs='+', help="kegg org code")
    parser_b = subparsers.add_parser('fq', help='input is fastq format')
    parser_b.add_argument('-b', choices='xyz', help='baz help')

    args = parser.parse_args()
    print("args after parse_args:", args, "type", type(args))
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
    params = Params()
    try:
        if argv is None:
            argv = sys.argv
        args = params.parse_options(argv)
        print(args)
    except Usage as err:
        print(sys.stderr, err.msg)
        print(sys.stderr, "for help use --help")
        return 2

if __name__ == "__main__":
    sys.exit(main())
