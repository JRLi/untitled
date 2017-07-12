#!/usr/bin/env python3.6
import os
import shutil


def main():
    anno_count, dir_count= 0, 0
    for dirPath, dirNames, fileNames in os.walk('.'):
        if dirPath.endswith('miRNA'):
            dir_count += 1
            print(dirPath)
            for f in fileNames:
                if f.endswith('.fasta'):
                    print('\t', f)
                    anno_count += 1
                    shutil.copy(os.path.join(dirPath, f),'default_fa/' + f.replace('fasta', 'fa'))
    print('dir:{}\tfa:{}'.format(dir_count, anno_count))

if __name__ == '__main__':
    main()