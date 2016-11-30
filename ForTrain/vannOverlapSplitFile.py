# -*- coding: UTF-8 -*-
"""Use for specific orchid annotation file: vannChart transcript to expression and annotation of refSeq results"""
import os

path = "D:/Project/orchid/orchid_sim90_len500_FDR01_report/vannOverlap/"    # vann output overlap file dir
path2 = "D:/Project/orchid/orchid_sim90_len500_FDR01_report/"
dirFiles = os.listdir(path)
print("dirFiles", dirFiles)

# dict 1: 7 area of 3 type of vannCharts, total 21 sub-area and transcripts
subFile2transListDict = dict()
for fileName in dirFiles:
    print("filename", fileName)
    inputFile = open(path+fileName, "r")
    subAreaName = ""
    for line in inputFile.readlines():
        line = line.replace("\"", "")
        print(line)
        if line.startswith("Lip"):
            if subAreaName != "":
                print(subAreaName, "end")
                subFile2transListDict[subAreaName] = transList
            subAreaName = line.replace("\n", "").replace("\"", "").replace("$", "_")
            transList = []
            print(subAreaName, "first step")
        else:
            stringField = line.replace("\"", "").replace("\n", "").split(sep=" ")
            print(stringField[0])
            transList.append(stringField[1])
    subFile2transListDict[subAreaName] = transList
    inputFile.close()

# dict 2: transcripts to expression and annotation of refSeq
trans2expressRefDict = dict()
expressRefFile = open(path2+"expressAndRefSeq.txt", "r")
for line in expressRefFile.readlines():
    if line.startswith("trans"):
        continue
    stringField = line.replace("\"", "").replace("\n", "").split(sep="\t", maxsplit=1)
    transID = stringField[0]
    expressRef = stringField[1]
    trans2expressRefDict[transID] = expressRef
expressRefFile.close()

# dict 3: transcripts to blast NR/NT result
'''
wait for blastx and blastn
'''
id2nt = open(path2+"orchid_scaffolds_no_extension_90_500.fa.nt.w28e-3_3_z_nt_name", "r")
id2ntDict = {}
lineCount = 0
for line in id2nt.readlines():
    lineCount += 1
    if lineCount == 1:
        continue
    stringField = line.replace("\n", "").split("\t")
    annotationList = id2ntDict.get(stringField[0], [])
    annotationString = stringField[1] + "\t" + stringField[2] + "\t" + stringField[7] + "\t" + stringField[8] + "\t"
    annotationList.append(annotationString)
    id2ntDict[stringField[0]] = annotationList
id2nt.close()

for outFileName in subFile2transListDict.keys():
    outputFile = open(path2+outFileName, "w")
    outputFile.write("transcript_ID\tlength\tLH3-Lip-TPM\tLH3-P-TPM\tLH6-Lip-TPM\tLH6-P-TPM\tLH-MF-Lip-TPM\t"
                     "LH-MF-P-TPM\tLH3-Lip-FPKM\tLH3-P-FPKM\tLH6-Lip-FPKM\tLH6-P-FPKM\tLH-MF-Lip-FPKM\tLH-MF-P-FPKM\t"
                     "Blast_to_RefSeq\tE-value\tBit Score\tblast_NT_1\tNT_annotation_1\tE-value_1\tBit-score_1\t"
                     "blast_NT_2\tNT_annotation_2\tE-value_2\tBit-score_2\tblast_NT_3\tNT_annotation_3\tE-value_3\t"
                     "Bit-score_3\n")
    for transcript in subFile2transListDict[outFileName]:
        outputFile.write(transcript+"\t"+trans2expressRefDict[transcript]+"\t")     # write annotation of transcript
        outList = id2ntDict.get(transcript, [])
        for line in outList:
            outputFile.write(line)
        outputFile.write("\n")
    outputFile.close()
