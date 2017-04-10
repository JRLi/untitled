#!/usr/bin/env python3.6
import os
import xml.etree.ElementTree as ET

path='D://Project/drs/clinical_test/'

def main():
    dir_count, all_count, nxml_count, gxml_count, txlsx_count, pass_count = 0, 0, 0, 0, 0, 0
    for dirPath, dirNames, fileNames in os.walk(path):
        if (dirPath == path) or (dirPath.endswith('logs')):
            continue
        dir_count += 1
        for f in fileNames:
            all_count += 1
            file_path = os.path.join(dirPath, f)
            print(file_path)
            if f.startswith('nationwide'):
                nxml_count+= 1
                tree = ET.parse(file_path)
                root = tree.getroot()
                print(root.tag)
                print(root.attrib)
            elif f.startswith('genome'):
                gxml_count += 1
            elif f.endswith('.xlsx'):
                txlsx_count += 1
            else:
                pass_count += 1
    print('dir: {}\nall: {}\nnxml: {}'.format(dir_count, all_count, nxml_count))
    print('gxml: {}\ntxlsx: {}\npass: {}'.format(gxml_count, txlsx_count, pass_count))

if __name__ == '__main__':
    main()