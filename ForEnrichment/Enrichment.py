# -*- coding = utf-8 -*-

import sys
from collections import defaultdict
from ForEnrichment import geneAnnot

pathOfPathway = '~/gene/'
pathOfNCBIGeneInfo = '~/gene/gene.May2015/'
orgID = '9606'
geneid2godict = defaultdict(lambda: set())
go2geneiddict = defaultdict(lambda: set())

if len(sys.argv) != 4:
    sys.exit("Usage: python3.5 Enrichment inputFilename pValue_th option")
inputfileName = sys.argv[1]
p_value_th = sys.argv[2]
option = sys.argv[3]


def enrich():
    geneAnnot.parse_gene_info(orgID)
    if option == "GO":
        go_analysis()
    if option == "pathway":
        print()


def go_analysis():
    parse_gene2go(orgID)
    type2genedict = parse_data(inputfileName)


def parse_gene2go(tax_id):
    genecount = 0
    with open(pathOfNCBIGeneInfo + 'gene2go.' + tax_id) as inputfile:
        for line in inputfile:
            items = line.replace("\n", "").replace("\"", "").split("\t")
            if len(items) < 2 or items[0] != tax_id:
                continue
            genecount += 1
            geneid = items[1]
            goid = items[2] + " " + items[5]
            goset = geneid2godict[geneid]
            goset.add(goid)
            geneidset = go2geneiddict[goid]
            geneidset.add(geneid)
    print("process:", pathOfNCBIGeneInfo + 'gene2go.' + tax_id)
    print("gene2go database, geneCount:", genecount)
    print("gene2go database, geneid2godict_size:", len(geneid2godict))
    print("gene2go database, go2geneiddict_size:", len(go2geneiddict))
    return go2geneiddict


def parse_data(filename):
    type2genedict = defaultdict(lambda: set())
    with open("./" + filename) as inputfile:
        for line in inputfile:
            items = line.replace("\n", "").replace("\"", "").split("\t")
            if len(items) < 2:
                continue
            tp = items[0]
            gene = items[1]
            sym = geneAnnot.getgenesym(gene)
            if sym == '-':
                gene = geneAnnot.getgeneid(gene)
            sym = geneAnnot.getgenesym(gene)
            if sym == '-':
                continue
            genes = type2genedict[tp]
            genes.add(sym)
        print("parse_data type2gene:", len(type2genedict))
    return type2genedict

enrich()







