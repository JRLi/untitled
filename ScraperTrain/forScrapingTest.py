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

htmlTitle = get_title(targetUrl)
if htmlTitle:
    print(htmlTitle)
else:
    print('Title could not be found')
