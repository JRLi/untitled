"""
usage: python3.5 decreaseBlastResult.py blastOutputFile(fileName, end with e-value) maxTarget(int)
Must be used before blastAnnotMin
"""
import sys
# control argv length, 3 or 4, if > 4 or < 3, print "usage: python .py blastFile int targetAnnotation(option)"
blastRawFile = str(sys.argv[1])
blastNum = int(sys.argv[2])
inputFile = open("./" + blastRawFile, "r")
outputFile = open("./" + blastRawFile + "_" + str(blastNum), "w")

queryID = ""
targetCheck = ""
blastDict = dict()
for line in inputFile.readlines():
    stringField = line.replace("\n", "").split(sep="\t", maxsplit=2)
    otherInf = stringField[1] + "\t" + stringField[2]
    if stringField[0] != queryID:
        if queryID != "":
            blastDict[queryID] = infList
        queryID = stringField[0]
        infList = [otherInf]
        targetCheck = stringField[1]
    else:
        if targetCheck == stringField[1]:
            continue
        infList.append(otherInf)
        targetCheck = stringField[1]
        '''
        if stringField[1] != targetCheck:
            infList.append(otherInf)
            targetCheck = stringField[1]
        else:
            continue
            '''
blastDict[queryID] = infList

inputFile.close()

forCount = 0
for outID in blastDict.keys():
    for infContext in blastDict[outID]:
        if forCount >= blastNum:
            forCount = 0
            break
        else:
            outputFile.write(outID+"\t"+infContext+"\n")
            forCount += 1

outputFile.close()
