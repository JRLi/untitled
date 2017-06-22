#!/usr/bin/env python3.6
import os

description_column = 2


def file_processor(file_path, m_c):
    with open(file_path) as inputFile:
        snv_count = 0; more2 = 0; more3 = 0; more5 = 0; more10 = 0; patients = 0
        for line in inputFile:
            lf = line.strip().split(',')
            if line.startswith(','):
                groups = ['-'.join(x.split('-')[:4]) for x in lf[m_c:]]
                patients = len(groups)
            else:
                snv_count += 1
                check = 0
                for x in lf[m_c:]:
                    if x != '0':
                        check += 1
                if check >= 2:
                    more2 += 1
                if check >= 3:
                    more3 += 1
                if check >= 5:
                    more5 += 1
                if check >= 10:
                    more10 += 1
    return patients, snv_count, more2, more3, more5, more10


def main():
    file_list = os.listdir('.')
    with open('status.txt', 'w') as out_status:
        out_status.write('snv_file\tpatients\ttotal_snv\tmore2\tmore2rate\tmore3\tmore3rate\t'
                         'more5\tmore5rate\tmore10\tmore10rate\n')
        for file in file_list:
            if file.startswith('snv_indel.extracted'):
                patient, snv, m2, m3, m5, m10 = file_processor(file, description_column)
                print(patient, snv, m2, m3, m5, m10)
                out_status.write('{}\t{}\t{}\t{}\t{:.3f}\t{}\t{:.4f}\t{}\t{:.5f}\t{}\t{:.8f}\n'.format(
                    file, patient, snv, m2, m2/snv, m3, m3/snv, m5, m5/snv, m10, m10/snv))

if __name__ == "__main__":
    main()