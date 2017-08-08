#!/usr/bin/env python
import os

path = 'D:/Project/circ_miRNA/barcode_files'
out_path = 'D:/Project/circ_miRNA/out_annotation'
def main():
    file_list = os.listdir(path)
    file_count = 0
    with open(os.path.join(out_path, 'all_mirna_ann.csv'), 'w') as out_f:
        for file in file_list:
            if file.startswith('T'):
                file_count += 1
                with open(os.path.join(path, file)) as in_f, open(os.path.join(out_path, 'lab_' + file), 'w') as lab_f:
                    line1 = next(in_f).rstrip()
                    n_count = 0
                    lab_f.write('cases,label\n')
                    if file_count == 1:
                        out_f.write(line1 + '\n')
                    for line in in_f:
                        lf = line.rstrip().split(',')
                        if lf[2] != '':
                            out_f.write(line.rstrip() + '\n')
                            check = lf[2].split('-')[3][0:2]
                            label = '1' if int(check) <= 9 else '0'
                            lab_f.write(','.join([lf[2], label]) + '\n')
                        else:
                            n_count += 1
                            lf[2] = 'unknown_case_' + str(n_count)
                            out_f.write(','.join(lf) + '\n')
                            lab_f.write(','.join([lf[2], '2']) + '\n')

                    if n_count != 0:
                        print('{}: {} file IDs have no case barcode'.format(file, n_count))
        print(file_count)

if __name__ == '__main__':
    main()
