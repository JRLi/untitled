#!/usr/bin/env python3.6
import os

def main():
    anno_count, dir_count= 0, 0
    for dirPath, dirNames, fileNames in os.walk('.'):
        if (dirPath in ['.', './miRNA_exp', './processed_coll']) or (dirPath.endswith('miRNA')):
            continue
        dir_count += 1
        print(dirPath, fileNames)
        for f in fileNames:
            if f.endswith('.collapse.fa'):
                anno_count += 1
                with open(os.path.join(dirPath, f)) as in_coll, open('processed_coll/' + f, 'w') as outFile:
                    for line in in_coll:
                        line = line.replace('-', '\t') if line.startswith('>') else line
                        outFile.write(line)
    print(dir_count, anno_count)

if __name__ == '__main__':
    main()