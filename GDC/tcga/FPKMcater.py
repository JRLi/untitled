#!/usr/bin/env python3.6
import os
import pandas as pd


def main():
    with open('./annotation.txt', 'w') as anno_file:
        df2 = pd.DataFrame()
        gz_count, anno_count, dir_count, all_count = 0, 0, 0, 0
        for dirPath, dirNames, fileNames in os.walk('.'):
            if (dirPath == '.') or (dirPath.endswith('logs')):
                continue
            dir_count += 1
            for f in fileNames:
                all_count += 1
                file_path = os.path.join(dirPath, f)
                fpath, uuid = os.path.split(dirPath)
                if f.startswith('annotation'):
                    anno_count += 1
                    with open(file_path) as inputFile:
                        for line in inputFile:
                            if line.startswith('id\t') and anno_count != 1:
                                continue
                            anno_file.write(line)
                elif f.endswith('FPKM.txt.gz'):
                    gz_count += 1
                    df1 = pd.read_table(file_path, index_col=0, header=None, names=[uuid])
                    df2[uuid] = df1[uuid]
                else:
                    pass
        print('df2.shape:', df2.shape)
        df2.to_csv('./rna_seq_FPKM.csv')
        print('FPKM.txt.gz: {}\nannotate: {}\ndir: {}\nall_process: {}\n'.format(gz_count, anno_count, dir_count, all_count))

if __name__ == '__main__':
    main()