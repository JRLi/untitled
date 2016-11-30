# -*- coding = utf-8 -*-
"""
usage: python3.5 ParseKeggGene.py hsaIDName.txt
Whatever filename you give the hsaIDName.txt, the format must be given in: ID \t Name \n.
No need to change anything unless the length of organism code is not 3 or the pathwayID contain another org code.
"""

import sys
import re

ID2PathwayFilename = str(sys.argv[1])
wwwFileCount = 0

with open("./" + ID2PathwayFilename) as inputID2Pathway, \
 open("./" + ID2PathwayFilename + "_ParseKeggGene", "w") as output1,\
 open("./" + ID2PathwayFilename + "_gene2KEGG", "w") as output2:
    for line in inputID2Pathway:
        stringField = line.replace("\"", "").replace("\n", "").split("\t")
        pathID = stringField[0]
        pathway = stringField[1]
        orgID = pathID[0:3]
        wwwFileName = 'www_bget?pathway:' + pathID

        with open("./" + wwwFileName) as wwwFile:
            wwwFileCount += 1
            # print("process:", wwwFileName)
            wwwLine = wwwFile.read().replace("\n", "")
            pattern = '\/dbget-bin\/www_bget\?' + orgID + ':\d+'
            match = re.findall(pattern, wwwLine)
            if match:
                for geneString in match:
                    geneID = geneString[geneString.rfind(":") + 1:]
                    output1.write(pathID + "\t" + pathway + "\t" + geneID + "\n")
                    output2.write(geneID + "\t" + pathID + " " + pathway + "\n")
            else:
                print(wwwFileName, "have no geneID")

print("Number of wwwFile processed:", wwwFileCount)
