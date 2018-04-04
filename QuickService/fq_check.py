#!/usr/bin/env python3.6
import sys


def main():
    with open(sys.argv[1]) as if1, open('conflict_read', 'w') as of1:
        lc = 0
        cf_dict = {}
        r_name, r_seq, r_qt = '', '', ''
        for line in if1:
            lc += 1
            if lc % 4 == 1:
                r_name = line.rstrip()
            if lc % 4 == 2:
                r_seq = line.rstrip()
            if lc % 4 == 0:
                r_qt = line.rstrip()
                if len(r_seq) != len(r_qt):
                    print('Line: {}'.format(lc))
                    print('{}: {} vs {}'.format(r_name, len(r_seq), len(r_qt)))
                    cf_dict[r_name + '\n' + r_seq ] = r_qt
        for k, v in cf_dict.items():
            of1.write('{}\n{}\n'.format(k, v))


if __name__ == 'main':
    main()