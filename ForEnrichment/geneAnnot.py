pathOfPathway = '~/gene/'
pathOfNCBIGeneInfo = '~/gene/gene.May2015/'
orgID = '9606'
geneid2symdict = {}
geneid2descdict = {}
sym2geneiddict = {}
geneid2aliasdict = {}
alias2symdict = {}


def parse_gene_info(tax_id):
    with open(pathOfNCBIGeneInfo + "gene_info." + tax_id) as geneInfoFile:
        print('parse_gene_info:', pathOfNCBIGeneInfo + "gene_info." + tax_id)
        gene_count = 0
        for line in geneInfoFile:
            if line[:line.find("\t") - 1] != tax_id:
                continue
            s_field = line.replace("\n", "").replace("\"", "").split("\t")
            if len(s_field) <= 9:
                print('length is less than 9:', line)
                continue
            gene_id = s_field[1]
            symbol = s_field[2]
            locus = s_field[3]
            synonyms = s_field[4]
            description = s_field[8]
            if len(gene_id) <= 1 or len(symbol) <= 1:
                continue
            gene_count += 1
            geneid2symdict[gene_id] = symbol
            geneid2descdict[gene_id] = description
            sym2geneiddict[symbol] = gene_id
            alias = set()
            alias.add(symbol)
            alias2symdict[symbol.upper()] = symbol
            alias.add(symbol.upper())
            items = synonyms.replace("|", " ").split(" ")
            if len(locus) >= 2:
                alias2symdict[locus] = symbol
                alias2symdict[locus.upper()] = symbol
                alias.add(locus)
                alias.add(locus.upper())
            for item in items:
                alias2symdict[item] = symbol
                alias2symdict[item.upper()] = symbol
                alias.add(item)
                alias.add(item.upper())
            if len(description) < 20 and description.endswith("p"):
                alias.add(description)
            geneid2aliasdict[gene_id] = alias
    print('geneid2symdict size:', len(geneid2symdict))
    print('sym2geneiddict_size:', len(sym2geneiddict))
    print('alias2symdict_size:', len(alias2symdict))
    print('geneid2descdict_size:', len(geneid2descdict))
    return sym2geneiddict


def getgenesym(geneid):
    return geneid2symdict.get(geneid, "-")


def getgeneid(sym):
    geneid = str(sym2geneiddict.get(sym, '-'))
    if geneid == '-':
        sym2 = str(alias2symdict.get(sym, '-'))
        geneid = str(sym2geneiddict.get(sym2, '-'))
    return geneid
print(len(geneid2symdict))