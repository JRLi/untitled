import os

file_path = 'D:/Project/drs/GDC/query_project/'
file_list = os.listdir(file_path)
print(len(file_list))

file_count = 0

with open(file_path + 'all_37_proj_fpkm.csv', 'w') as out_file:
    for file in file_list:
        file_count += 1
        with open(file_path + file) as in_file:
            for line in in_file:
                if line.startswith('file_id') and file_count != 1:
                    continue
                out_file.write(line)