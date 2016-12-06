from urllib.request import urlopen
from urllib.error import HTTPError
from bs4 import BeautifulSoup

targetUrl = "http://www.pythonscraping.com/pages/page1.html"


def get_title(url):
    try:
        html = urlopen(url)
    except HTTPError as e:
        print(e)
        return None
    try:
        bsObj = BeautifulSoup(html.read(), "html.parser")
        title = bsObj.body.h1
    except AttributeError as e:
        print(e)
        return None
    return title
#
# htmlTitle = get_title(targetUrl)
# if htmlTitle:
#     print(htmlTitle)
# else:
#     print('Title could not be found')


def get_path_list(url):
    print("to get path list")
    try:
        html = urlopen(url)
        path_data = html.read().decode('utf-8')
        print("get list done")
    except HTTPError as e:
        print(e)
        return None
    return path_data
v_url = 'http://rest.kegg.jp/list/hsa'
gene_annotation_dict = {}
entry_data = get_path_list(v_url)
line_list = str(entry_data).replace("\"", "").replace('hsa' + ":", "").replace(" | ", "\t").split("\n")
for line in line_list:
    if len(line) > 3:
        items = line.split("\t", 1)
        gene_annotation_dict[items[0].strip()] = items[1].strip()
print(len(gene_annotation_dict))
print(gene_annotation_dict.get('100287010', ""))
