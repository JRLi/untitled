# -*- coding: UTF-8 -*-
"""
usage: python3.5 getBlastAnnotation.py blastDecreasedFile(fileName, end with e-value) TargetAnnotation(all line >)
Must be used after decreaseBlastResult.py
Must be used in linux server with large ram
If the target DB is a small single sequence file, e.g. refSeq file, use the java BlastAnnot.java instead.
"""
# target seq format: >ID(HQ832783.1) description
# blast file format: queryID \t target(gi|329755407|gb|HQ832783.1|) \t other information
# goal: queryID \t target(HQ832783.1 \t description) \t other information
# inf: queryStart queryEnd targetStart targetEnd Expect_value Bit_score gap_open identical_matches Alignment_length
import sys

# check, if use argv, may use the path rather then name
blastDecreasedFile = str(sys.argv[1])           # end with a _number
blastTargetAnnotationFile = str(sys.argv[2])    # no sequence, just >id description per line
inputFileBlast = open("./" + blastDecreasedFile, "r")
inputFileTarget = open("./" + blastTargetAnnotationFile, "r")
outputFile = open("./" + blastDecreasedFile + "_" + blastTargetAnnotationFile, "w")

# Construct dictionary between targetID and annotation description
targetID2annotationDict = {}
modID2annotationDict = {}
print("Start to construct the annotation dictionary")
lineCount = 0
for line in inputFileTarget.readlines():
    lineCount += 1
    stringField = line.replace("\n", "").split(sep=" ", maxsplit=1)
    targetID = stringField[0][1:]
    modID = stringField[0][1:-1]
    annotationInf = stringField[1]
    targetID2annotationDict[targetID] = annotationInf
    modID2annotationDict[modID] = annotationInf
    if lineCount % 100000 == 0:
        print(lineCount)
inputFileTarget.close()

print("The number of lines of the annotation file:", lineCount)
print("The number of dictionary items:", len(targetID2annotationDict), len(modID2annotationDict))
outputFile.write("transcript_ID\tblastn_target\tblastn_description\tStart_of_query\tEnd_of_query\tStart_of_target\t"
                 "End_of_target\tE-value\tBit-score\tidentical_matches\tAlignment_length\n")

lineCount = 0
for line in inputFileBlast.readlines():
    lineCount += 1
    stringField = line.replace("\n", "").split(sep="\t", maxsplit=2)
    giResult = stringField[1]
    if giResult.rfind("|") == (len(giResult) - 1):
        blastTarget = giResult[giResult.rfind("|", 0, -1) + 1:giResult.rfind("|")]
    else:
        blastTarget = giResult[giResult.rfind("|", 0, giResult.rfind("|")) + 1:giResult.rfind("|")]
    blastDescription = targetID2annotationDict.get(blastTarget, "-")
    if blastDescription == "-":
        blastTargetNew = blastTarget[:-1]
        blastDescription = modID2annotationDict.get(blastTargetNew, "-")
    outputFile.write(stringField[0] + "\t" + blastTarget + "\t" + blastDescription + "\t" + stringField[2] + "\n")
inputFileBlast.close()

print("The number of lines of the blast file:", lineCount)
outputFile.close()