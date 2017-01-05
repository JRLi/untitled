#!/usr/bin/env python3.5
from urllib.request import urlopen
from urllib.error import HTTPError
from bs4 import BeautifulSoup
from collections import OrderedDict
import itertools
import re
import os
import sys
import argparse

use_message = '''
    To get NCBI GEO information according GSM.
    Need Python3 and BeautifulSoup4.
    GSM example: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1000001
'''

gsmUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM"
gsmStart = 1000001  # use argparse to define this number
targetUrl = "http://www.genome.jp/dbget-bin/www_bget?pathway:"
output_dir = "./gene_pathway_out/"


class Gsm():
    def __init__(self, gsm, tax, char, title, desc, source, molecule):
        self.gsm = gsm
        self.tax = tax
        self.char = char
        self.title = title
        self.desc = desc
        self.source = source
        self.molecule = molecule


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


def get_html(url):
    print("to get gsm:", url)
    try:
        html = urlopen(url)
        html_data = html.read().decode('utf-8')
    except HTTPError as e:
        print(e)
        return None
    return html_data


def get_gsm_information(gsm_html_data, gsm_id):
    if gsm_html_data is None:
        print('There is no results of GSM_id:', 'GSM' + gsm_id)
        sys.exit()
    else:
        tag1, tag2 = "div", "a"
        regex_1 = '\/dbget-bin\/www_bget\?' + gsm_id + ':\w*\d+'
        regex_2 = ''
        print("Process GSM data:", 'GSM' + gsm_id)
        try:
            bs_object = BeautifulSoup(gsm_html_data)
            tax_id = bs_object.find(tag1).find_all(tag2, href=re.compile(regex_1))
            characteristics = bs_object.find(tag1).find_all(tag2, href=re.compile(regex_2))
            title = ''
            description = ''
            source_name = ''
            molecule = ''
        except AttributeError as e:
            print(e)
            return None
    return Gsm(gsm_id, tax_id, characteristics, title, description, source_name, molecule)

if __name__ == "__main__":
    # for i in range(sys.maxsize)
    for i in itertools.count(start=gsmStart):
        process_url = gsmUrl + str(i)
        gsm_html = get_html(process_url)
        gsm = get_gsm_information(gsm_html, i)
