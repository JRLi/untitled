# -*- coding: UTF-8 -*-
import os

path = "D:/Project/PBMC/finalCLS/total/"
dirFiles = os.listdir(path)
print(len(dirFiles))
print(dirFiles)
fileCount = 0

with open("E:/StringTemp/CLSs", "w") as outputFile:
    outputFile.write("GSE,GPL,Type,GSM,immuneCell,CLS\n")
    for fileName in dirFiles:
        fileField = fileName.replace(".txt", "").split("_")
        GSE = fileField[0]
        GPL = fileField[1]
        profileType = fileField[3]
        with open(path + fileName) as inputFile:
            fileCount += 1
            for line in inputFile:
                if line.startswith("\t"):
                    immuneList = line.replace("\n", "").split("\t")[1:-1]
                    print(fileCount, GSE, GPL, profileType, len(immuneList), immuneList)
                    CLSCount = 1
                else:
                    stringField = line.replace("\n", "").split("\t")
                    GSM = stringField[0]
                    CLSList = stringField[1:]
                    if CLSCount == 1:
                        print(len(CLSList))
                    CLSCount = 0
                    for i in range(len(CLSList)):
                        outputFile.write(GSE + "," + GPL + "," + profileType + "," + GSM + "," + immuneList[i] + "," +
                                         CLSList[i] + "\n")


