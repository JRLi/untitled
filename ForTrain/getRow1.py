# -*- coding: UTF-8 -*-

import os
path = "D://Project/PBMC/ImmuneProfile/"    # attention, the output file will present in the same dir
dirFiles = os.listdir(path)
print("dirFiles:", dirFiles)

geneIDSet = set()
for fileName in dirFiles:
    print("inputFile:", fileName)
    inputFile = open(path + fileName, "r")
    lineCount = 0
    for line in inputFile.readlines():
        lineCount += 1
        if lineCount == 1:
            continue
        lineList = line.replace("\n", "").replace("\"", "").split(sep="\t")
        geneIDSet.add(lineList[0])
    print("gene number:", lineCount - 1, len(geneIDSet))
    inputFile.close()

outputFile = open(path + "all_geneID", "w")
for geneID in geneIDSet:
    if geneID != "":
        outputFile.write(geneID+"\n")
outputFile.close()
