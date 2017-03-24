#!/usr/bin/env python3.5
import os
import sys
import numpy as np
import pandas as pd
import argparse
use_message = '''
    To get geneID kegg pathway annotation.
    Need Python3 and BeautifulSoup4.
    organism_code[e.g. hsa eco ecj]: http://rest.kegg.jp/list/organism
'''


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


def args_parse():
    parser = argparse.ArgumentParser(description=use_message)
    parser.add_argument('profile', nargs='+', help="kegg org code")
    parser.add_argument("-v", "--verbose", action="store_true", help="output verbose results")
    args = parser.parse_args()
    return args


def main(argv=None):
    try:
        if argv is None:
            argv = args_parse()
            file_list = argv.profile
            for profile in file_list:
                pass
    except Usage as err:
        print(sys.stderr, err.msg)
        print(sys.stderr, "Terminated, for help use -h or --help")
        return 2

if __name__ == '__main__':
    sys.exit(main())