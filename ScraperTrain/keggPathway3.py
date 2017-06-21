#!/usr/bin/env python3.6
"""
Created by JRLi on 2016/12/03 for python learning
"""

from urllib.request import urlopen
from urllib.error import HTTPError
from bs4 import BeautifulSoup
from collections import OrderedDict
import re
import os
import sys
import argparse

use_message = '''
    To get geneID kegg pathway annotation.
    Need Python3 and BeautifulSoup4.
    organism_code[e.g. hsa eco ecj]: http://rest.kegg.jp/list/organism
'''
pathAPIUrl = "http://rest.kegg.jp/list/pathway/"
geneAPIUrl = "http://rest.kegg.jp/list/"
targetUrl = "http://www.genome.jp/dbget-bin/www_bget?pathway:"
output_dir = "./gene_pathway_out/"
tag1 = "div"
tag2 = "a"


def args_parse():
    parser = argparse.ArgumentParser(description=use_message)
    parser.add_argument('organism', nargs='+', help="kegg org code")
    parser.add_argument("-v", "--verbose", action="store_true", help="output verbose results")
    args = parser.parse_args()
    return args


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
        return None
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


def get_verbose_annotation_dict(v_url, org):
    gene_annotation_dict = {}
    entry_data = get_path_list(v_url)
    line_list = str(entry_data).replace("\"", "").replace(org + ":", "").replace(" | ", "\t").split("\n")
    for line in line_list:
        if len(line) > 3:
            items = line.split("\t", 1)
            gene_annotation_dict[items[0].strip()] = items[1].strip()
    return gene_annotation_dict


def get_gene_path_list(pathwaydict, regex, verbose, vdict):
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
                    if verbose:
                        something = geneID + "\t" + vdict.get(geneID) + "\t" + pathid + " " + pathwaydict[pathid]
                    else:
                        something = geneID + "\t" + pathid + " " + pathwaydict[pathid]
                    gene_path_list.append(something)
            except AttributeError as e:
                print(e)
                return None
    return gene_path_list


def main(argv=None):
    if argv is None:
        argv = args_parse()
        print('Number of organism:', len(argv.organism), '\n', 'Detail:', argv.organism)
        print('Verbose mode' if argv.verbose else 'Default quiet mode')
        suffix = '.v' if argv.verbose else ''

    report_list = []
    for org_code in argv.organism:
        pattern = '\/dbget-bin\/www_bget\?' + org_code + ':\w*\d+'
        path_url = pathAPIUrl + org_code
        prepare_output_dir()
        pathway_data = get_path_list(path_url)
        pathway_name_dict = get_pathway_dict(pathway_data, org_code)
        if pathway_name_dict is None:
            report_list.append(org_code)
            continue
        print(len(pathway_name_dict))

        if argv.verbose:
            gene_list_url = geneAPIUrl + org_code
            gene_verbose_dict = get_verbose_annotation_dict(gene_list_url, org_code)
            print(org_code, "have", len(gene_verbose_dict), "entire.")
        else:
            gene_verbose_dict = {}

        gene_pathway_list = get_gene_path_list(pathway_name_dict, pattern, argv.verbose, gene_verbose_dict)

        if gene_pathway_list is None:
            print("There are some errors when getting gene_pathway_list")
            sys.exit()
        else:
            print("gene_pathway_list size:", len(gene_pathway_list))
            with open(output_dir + "gene2KEGG." + org_code + suffix, "w") as outputFile:
                for line in gene_pathway_list:
                    outputFile.write(line + "\n")
                print("Output file:", output_dir + "gene2KEGG." + org_code + suffix)
    print("No path url:", report_list)
    
if __name__ == "__main__":
    sys.exit(main())
