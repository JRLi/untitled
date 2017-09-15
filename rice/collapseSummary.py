import os
collapse_d = 'unmatch_AfterBlastmiRNA_formirplant'

def main():
    file_list = os.listdir(collapse_d)
    with open('collapse_status_rest', 'w') as out_f:
        out_f.write('Rice_Cultivar\tread_type\ttotal_read_count\n')
        for file in file_list:
            read_c, rc_total = 0, 0
            with open(os.path.join(collapse_d, file)) as in_f:
                for line in in_f:
                    if line.startswith('>'):
                        read_c += 1
                        lf = line.rstrip().split('-')
                        rc_total += int(lf[1])
                out_f.write('{}\t{}\t{}\n'.format(file.replace('_unmatched_collapse_m.fa', ''), read_c, rc_total))

if __name__ == '__main__':
    main()