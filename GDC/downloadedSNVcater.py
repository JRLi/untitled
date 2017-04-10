#!/usr/bin/env python3.6
import os
import gzip


def main():
    with open('./annotation.txt', 'w') as anno_file, open('./snv.maf', 'w') as snv_file:
        gz_count, anno_count, dir_count, all_count = 0, 0, 0, 0
        anno_file.write('id\tentity_type\tentity_id\tcategory\tclassification\tstatus\tnotes\n')
        for dirPath, dirNames, fileNames in os.walk('.'):
            if (dirPath == '.') or (dirPath.endswith('logs')):
                continue
            dir_count += 1
            for f in fileNames:
                all_count += 1
                file_path = os.path.join(dirPath, f)
                if f.startswith('annotation'):
                    anno_count += 1
                    with open(file_path, 'r') as inputFile:
                        for line in inputFile:
                            if line.startswith('id\t'):
                                continue
                            lf = line.replace('\"', '').replace('\r\n','').replace('\n', '').split('\t')
                            anno_file.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                                lf[0], lf[2], lf[3], lf[4], lf[5], lf[7], lf[8]))
                elif f.endswith('maf.gz'):
                    gz_count += 1
                    with gzip.open(file_path, 'rt') as inputFile:
                        for line in inputFile:
                            if (line.startswith('Hugo_Symbol') or line.startswith('#ver')) and gz_count != 1:
                                continue
                            snv_file.write(line)
                else:
                    pass
        print('maf.gz: {}\nannotate: {}\ndir: {}\nall_process: {}\n'.format(gz_count, anno_count, dir_count, all_count))

if __name__ == '__main__':
    main()