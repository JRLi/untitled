#!/usr/bin/env python3.6
import os, sys
import argparse


def args_parse():
    parser = argparse.ArgumentParser(description='for remove major zero row')
    parser.add_argument('-t', '--threshold', type=int, default=2, help="threshold of non-zero element in a row")
    parser.add_argument('file', nargs='?', help="file to process")
    args = parser.parse_args()
    return args


def prepare_output_dir(output_dir):
    print("prepare output dir")
    if os.path.exists(output_dir):
        pass
    else:
        os.mkdir(output_dir)


def file_processor(file_path, threshold, outdir):
    with open(file_path) as inputFile, open(outdir + '/' + file_path + '_' + str(threshold) + '.csv', 'w') as output:
        lineCount = 0
        for line in inputFile:
            lineCount += 1
            lf = line.strip().split(',')
            if lineCount == 1:
                groups = ['-'.join(x.split('-')[:4]) for x in lf[4:]]
                output.write(',' + ','.join(groups) + '\n')
            else:
                check = 0
                lf[1] = 'unknown' if lf[1] == '0' else lf[1]
                mut = ':'.join(lf[:4])
                for x in lf[4:]:
                    if x != '0':
                        check += 1
                if check >= threshold:
                    output.write(mut + ',' + ','.join(lf[4:]) + '\n')


def main(argv=None):
    if argv is None:
        argv = args_parse()
        print(argv)
        out_dir = 'contain' + str(argv.threshold)
        prepare_output_dir(out_dir)
        file_processor(argv.file, argv.threshold, out_dir)

if __name__ == "__main__":
    sys.exit(main())