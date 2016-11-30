from urllib.request import urlopen
from urllib.error import HTTPError
targetUrl = "http://rest.kegg.jp/list/pathway/"
orgCode = 'hsa'


def get_title(url):
    try:
        html = urlopen(url)
    except HTTPError as e:
        print(e)
        return None
    return html.read()

pathway_Name = get_title(targetUrl + orgCode)
if pathway_Name is None:
    print('There is no results of org_code:', orgCode)
else:
    print(get_title(targetUrl + orgCode))
