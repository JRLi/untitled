from collections import defaultdict

frequency = 50

def main():
    mut2cell_dict = defaultdict(list)
    cell_set = set()
    with open('D://Project/drs/gdsc/release_6.0/WES_variants.txt', 'r') as input_file, \
            open('D://Project/drs/gdsc/release_6.0/gdsc_WES_df_' + str(frequency) + '.txt', 'w') as output_file:
        line_count = 0
        for line in input_file:
            line_count += 1
            sf = line.split('\t')
            mut, cell_id = sf[1], sf[0]
            mut2cell_dict[mut].append(cell_id)
            cell_set.add(cell_id)
            check = line_count % 5000
            if check == 0:
                print(line_count)

        mut_list = list(mut2cell_dict)
        cell_list = list(cell_set)
        output_file.write('\t' + '\t'.join(cell_list) + '\n')

        for mut_id in mut_list:
            cell_check_list = mut2cell_dict[mut_id]
            if len(set(cell_check_list)) < frequency:
                continue
            v_list = [mut_id]
            for cell in cell_list:
                if cell in cell_check_list:
                    v_list.append('1')
                else:
                    v_list.append('0')
            output_file.write('\t'.join(v_list) + '\n')

if __name__ == '__main__':
    main()