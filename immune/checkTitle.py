#!/usr/bin/env python3.6
import os


def main():
    tsv_count, dir_count, title_set = 0, 0, set()
    for dirPath, dirNames, fileNames in os.walk('.'):
        if dirPath.startswith('./D00'):
            dir_count += 1
            tmp_count = 0
            for f in fileNames:
                if f.endswith('.tsv'):
                    title_check = len(title_set)
                    with open(os.path.join(dirPath, f)) as in_tsv:
                        tsv_count += 1
                        tmp_count += 1
                        f_line = in_tsv.readline().rstrip()
                        title_set.add(f_line)
                        if len(title_set) > title_check and len(title_set) > 1:
                            print('check:', dirPath, f)
            print(dirPath, tmp_count)
    print('dir:{}\ttsv:{}'.format(dir_count, tsv_count))
    print(len(title_set))
    for i in title_set:
        print(i)

if __name__ == '__main__':
    main()