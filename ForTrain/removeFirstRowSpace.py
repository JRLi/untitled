# -*- coding: UTF-8 -*-
# usage: python3.5 removeFirstRowSpace.py GSEFile.txt
import sys

inputFileName = str(sys.argv[1])
inputFile = open("./" + inputFileName, "r")
outputFile = open("./noSpaceGSE/" + inputFileName, "w")
lineCount = 0

for line in inputFile.readlines():
    lineCount += 1
    if lineCount == 1:
        outputFile.write(line)
    else:
        stringField = line.replace("\n", "").split(sep="\t", maxsplit=1)
        if stringField[0].startswith("Septin "):
            stringField[0] = stringField[0].replace("Septin ", "Sept")
        elif " " in stringField[0]:
            stringField[0] = stringField[0].replace(" ", "")
        outputFile.write(stringField[0] + "\t" + stringField[1] + "\n")

print("File:", inputFileName, "line_count:", lineCount)
inputFile.close()
outputFile.close()

