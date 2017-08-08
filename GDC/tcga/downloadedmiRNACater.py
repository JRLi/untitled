#!/usr/bin/env python3.6
import os
import pandas as pd


def main():
    with open('./annotation.txt', 'w') as anno_file:
        df_rd = pd.DataFrame()
        df_nr = pd.DataFrame()
        df_cr = pd.DataFrame()
        mir_count, anno_count, dir_count, all_count = 0, 0, 0, 0
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
                elif f.endswith('quantification.txt'):
                    mir_count += 1
                    df1 = pd.read_table(file_path, index_col=0)
                    df_rd[uuid] = df1['read_count']
                    df_nr[uuid] = df1['reads_per_million_miRNA_mapped']
                    df_cr[uuid] = df1['cross-mapped']
                else:
                    pass
        print('rd.shape: {}\tnr_shape: {}\tcr_shape: {}'.format(df_rd.shape, df_nr.shape, df_cr.shape))
        df_rd.to_csv('./mirna_read_count.csv')
        df_nr.to_csv('./mirna_rpm.csv')
        df_cr.to_csv('./mirna_cross_mapped.csv')
        print('mir: {}\nannotate: {}\ndir: {}\nall_process: {}\n'.format(mir_count, anno_count, dir_count, all_count))

if __name__ == '__main__':
    main()