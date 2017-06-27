# -*- coding = UTF-8 -*-
import os
overlapPath = "D:/Project/orchid/orchid_sim90_len500_FDR01_report/vannOverlap/"
faPath = "D:/Project/orchid/orchid_sim90_len500_FDR01_report/"
dirFiles = os.listdir(overlapPath)
print(dirFiles)

orchidID2SeqDict = {}
sequence, orchidID = "", ""
IDCount = 0
with open(faPath + "orchid_scaffolds_no_extension_90_500.fa") as inputFileOrchid:
    for line in inputFileOrchid:
        if line.startswith(">"):
            if orchidID != "":
                orchidID2SeqDict[orchidID] = sequence
            IDCount += 1
            orchidID = line.replace("\n", "")[1:]
            sequence = ""                       # VERY IMPORTANT, CRUCIAL STEP, RESET THE SEQ
        else:
            sequence += line.replace("\n", "")
    orchidID2SeqDict[orchidID] = sequence

print("orchid file transcripts:", IDCount)

orchidOverlapSet = set()
for fileName in dirFiles:
    with open(overlapPath + fileName) as inputFileOverlap:
        for line in inputFileOverlap:
            line = line.replace("\"", "").replace("\n", "")
            if not line.startswith("Lip"):
                stringField = line.split(sep=" ")
                orchidOverlapSet.add(stringField[1])

print("overlapSet_size:", len(orchidOverlapSet), " faDict_size:", len(orchidID2SeqDict))

with open(faPath + "overlapOrchid500.fa", "w") as outputFile:
    for oID in orchidOverlapSet:
        outputFile.write(">" + oID + "\n" + orchidID2SeqDict[oID] + "\n")
