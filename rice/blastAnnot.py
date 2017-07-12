#!/usr/bin/env python3
from collections import defaultdict
import sys
path = 'D:/Project/Rice/osa_mature.fa'


def id2seq_dict(path):
    with open(path) as in_fa:
        id, seq, id2seq, seq2id = '', [], {}, defaultdict(set)
        for line in in_fa:
            if line.startswith('>'):
                if id != '':
                    id2seq[id] = ''.join(seq)
                    seq2id[''.join(seq)].add(id)
                    seq = []
                id = line.rstrip().split(' ')[0].lstrip('>')
            else:
                seq.append(line.rstrip())
        return id2seq, seq2id


def main(argv=None):
    if argv is None:
        id2seq, seq2idSet = id2seq_dict(path)
        print(len(id2seq), len(seq2idSet))
        for k, v in seq2idSet.items():
            if len(v) > 1:
                print(k, v)

if __name__ == '__main__':
        sys.exit(main())
