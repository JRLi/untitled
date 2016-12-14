import os
from urllib.request import urlopen
from urllib.error import HTTPError
import csv

geneAPIUrl = "http://rest.kegg.jp/list/"
path='E:/StringTemp/new/'


def get_urlEntry(url):
    try:
        html = urlopen(url)
        entry_data = html.read().decode('utf-8')
        return entry_data
    except HTTPError as e:
        print(e)
        return None


def annotation_dict(url, org):
    entry_dict = {}
    entry_data = get_urlEntry(url)
    line_list = str(entry_data).replace("\"", "").replace(org + ":", "").replace(" | ", "\t").split("\n")
    for line in line_list:
        if len(line) > 3:
            items = line.split('\t')
            symbol = items[2][items[2].find(')') + 2:items[2].find(';')]
            entry_dict[symbol.strip()] = items[0].strip()
            entry_dict[symbol.strip().lower()] = items[0].strip()   # special for this module
    return entry_dict

fileList = os.listdir(path)
print(fileList)

for fileName in fileList:
    print('Now is the file: {} to process.'.format(fileName))
    orgID = 'ecj' if fileName.find('3110') > 0 else 'eco'
    prefix = 'K12_W3110' if fileName.find('3110') > 0 else 'K12_MG1655'
    target_url = geneAPIUrl + orgID
    symbol_to_id_dict = annotation_dict(target_url, orgID)
    line_count = 0
    with open(path + fileName, 'r') as inputFile,\
    open(path + prefix + orgID + fileName[fileName.find('_with'):fileName.find('_trans')] + '_1.csv', 'w', newline='') as out1,\
    open(path + prefix + orgID + fileName[fileName.find('_with'):fileName.find('_trans')] + '_2.csv', 'w', newline='') as out2:
        csv1_writer = csv.writer(out1)
        csv2_writer = csv.writer(out2)
        for line in inputFile:
            line_count += 1
            items = line.replace('\n', '').replace('\"', '').split('\t')
            if line_count == 1:
                c_normal, c_1, c_2 = items[2][5:], items[3][5:], items[4][5:]
                continue
            if len(items[1]) <= 2:
                continue
            k_id = symbol_to_id_dict.get(items[1],'-')
            e_normal = float(items[2]) if float(items[2]) >= 0.1 else 0.1
            e_1 = float(items[3]) if float(items[3]) >= 0.1 else 0.1
            e_2 = float(items[4]) if float(items[4]) >= 0.1 else 0.1
            f_1n = e_1/e_normal
            f_2n = e_2/e_normal
            csv1_writer.writerows([[c_1, c_normal, k_id, items[1],f_1n]])
            csv2_writer.writerows([[c_2, c_normal, k_id, items[1],f_2n]])
