# encoding: utf-8
"""
Created by JRLi on 2016/11/30 for python learning
"""

from urllib.request import urlopen
from urllib.error import HTTPError
from bs4 import BeautifulSoup
from collections import OrderedDict
import re
import os
import sys

use_message = '''
Need Python3 and BeautifulSoup
Usage:
    python3.5 keggPathway.py organism_code[e.g. hsa eco ecj]
    organism_code: http://rest.kegg.jp/list/organism
'''
pathAPIUrl = "http://rest.kegg.jp/list/pathway/"
targetUrl = "http://www.genome.jp/dbget-bin/www_bget?pathway:"
output_dir = "./gene_pathway_out/"
tag1 = "div"
tag2 = "a"


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


def prepare_output_dir():
    print("prepare output dir")
    if os.path.exists(output_dir):
        pass
    else:
        os.mkdir(output_dir)


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


def get_pathway_dict(pathway, org):
    if pathway is None:
        print('There is no results of org_code:', org)
        sys.exit()
    else:
        print("get pathway dict")
        pathway_dict = OrderedDict()
        line_list = str(pathway).replace("\"", "").replace("path:", "").split("\n")
        # print(line_list)
        for line in line_list:
            if len(line) > 3:
                items = line.split("\t")
                pathway_dict[items[0].strip()] = items[1].strip()
        print("get pathway dict done")
        return pathway_dict


def get_gene_path_list(pathwaydict, regex):
    if len(pathwaydict.keys()) == 0:
        print("There are some errors cause pathwaydict is empty")
        sys.exit()
    else:
        print("Numbers ot pathway:", len(pathwaydict.keys()))
        gene_path_list = []
        for pathid in pathwaydict.keys():
            print(pathid)
            try:
                html = urlopen(targetUrl + pathid)
            except HTTPError as e:
                print(e)
                return None
            try:
                bsObj = BeautifulSoup(html)
                target = bsObj.find(tag1).find_all(tag2, href=re.compile(regex))
                for line in target:
                    print(line)
                    geneID = line.string
                    something = geneID + "\t" + pathid + " " + pathwaydict[pathid]
                    gene_path_list.append(something)
            except AttributeError as e:
                print(e)
                return None
    return gene_path_list


def main(argv=None):
    if argv is None:
        argv = sys.argv
    if len(argv) != 2:
        raise Usage(use_message)
    else:
        orgCode = sys.argv[1]
        if orgCode in ('-h', '-help', '--help'):
            print(use_message)
            exit()
        pattern = '\/dbget-bin\/www_bget\?' + orgCode + ':\w*\d+'
        pathUrl = pathAPIUrl + orgCode
        prepare_output_dir()
        pathway_data = get_path_list(pathUrl)
        pathway_name_dict = get_pathway_dict(pathway_data, orgCode)
        print(len(pathway_name_dict))
        geneid_pathway_list = get_gene_path_list(pathway_name_dict, pattern)
        if geneid_pathway_list is None:
            print("There are some errors when getting geneid_pathway_list")
            sys.exit()
        else:
            print("geneid_pathway_list size:", len(geneid_pathway_list))
        with open(output_dir + "gene2KEGG." + orgCode, "w") as outputFile:
            for line in geneid_pathway_list:
                outputFile.write(line + "\n")
        print("Output file:", output_dir + "/gene2KEGG." + orgCode)

if __name__ == "__main__":
    sys.exit(main())


