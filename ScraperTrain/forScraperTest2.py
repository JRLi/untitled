from urllib.request import urlopen
from urllib.error import HTTPError
from bs4 import BeautifulSoup
import re

targetUrl = "http://www.genome.jp/dbget-bin/www_bget?pathway:hsa00010"
orgID = 'hsa'
pattern = '\/dbget-bin\/www_bget\?' + orgID + ':\d+'
def get_title(url):
    try:
        html = urlopen(url)
    except HTTPError as e:
        print(e)
        return None
    try:
        bsObj = BeautifulSoup(html)
        title = bsObj.find("div").find_all("a", href=re.compile(pattern))
    except AttributeError as e:
        print(e)
        return None
    return title

outFile = open("D:/Project/orchid/orchid_sim90_len500_FDR01_report/hsa00010_gene", "w")
htmlTitle = get_title(targetUrl)
if htmlTitle:
    print(type(htmlTitle))
    for line in htmlTitle:
        # geneID = line[line.rfind(">") + 1:line.rfind("<")]
        geneID = line.string
        outFile.write(geneID + "\n")
        print(type(line.string), line.string)
else:
    print('Title could not be found')
