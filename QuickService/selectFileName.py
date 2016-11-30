# -*- coding: UTF-8 -*-

import sys

typeFile = str(sys.argv[1])
inputFile = open("./" + typeFile, "r")
outputFileOne = open("./" + "ch1List", "w")
outputFileTwo = open("./" + "ch2List", "w")

for line in inputFile.readlines():
    stringField = line.replace("\n", "").split(" ")
    fileName = stringField[1].rstrip(",")
    print(fileName, stringField[-1])
    if stringField[-1] == 'one.channel':
        outputFileOne.write(fileName + " ")
    elif stringField[-1] == 'two.channel':
        outputFileTwo.write(fileName + " ")

inputFile.close()
outputFileOne.close()
outputFileTwo.close()
