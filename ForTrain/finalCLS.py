# -*- coding: UTF-8 -*-
""" To calculate the final CLS with CLSup/dn and output all results """
import os

hematopoieticType = ['Mouse', 'H1_GSE22886', 'H2_GSE24759']     # subDir name, must change if needed
fileCount = 0
for cellType in hematopoieticType:
    path = "D:/Project/PBMC/CLSupAnddn/" + cellType     # input PBMC CLSup/dn output file dir
    path2 = "D:/Project/PBMC/finalCLS/" + cellType      # output final CLS dir
    dirFiles = os.listdir(path)
    print("dirFiles:", dirFiles)
    for fileName in dirFiles:
        print("inputFile:", fileName)
        inputFile = open(path + "/" + fileName, "r")    # the file must in subDir
        outputFile = open(path2 + "/" + fileName + "_final", "w")   # must create subDir in the final CLS dir
        fileCount += 1
        for line in inputFile.readlines():
            if not line.startswith("GSM"):  # sample must come from GEO, if not, must change this area
                stringField = line.replace("_up.ES", "").replace("_dn.ES", "").split("\t")
                cellLength = len(stringField)//2
                print("up and down element:", len(stringField), "first cell:", stringField[0],
                      ", last cell:", stringField[len(stringField)//2 - 1])
                print("The first cell index with down regulator:", line.find(stringField[0], 1))
                outputFile.write("\t"+line[0:line.find(stringField[0], 1)].replace("_up.ES", "")+"\n")  # all cells ID
            else:
                pCLSField = line.split("\t")
                outputFile.write(pCLSField[0])  # start of per line with sample ID
                for j in range(1, cellLength + 1):
                    CLS = float(pCLSField[j]) - float(pCLSField[j + cellLength])
                    outputFile.write("\t"+str(CLS))
                outputFile.write("\n")  # end of per line
        inputFile.close()
        outputFile.close()
print("Total files processed:", fileCount)    # final check the number of file processed
