filePath = 'E:/StringTemp/GDS/'
fileNameIn = 'GDS3876.soft'
fileNameOut = 'GDS3876.matrix'

with open(filePath + fileNameIn, 'r') as raw_file, open(filePath + fileNameOut, 'w') as mod_file:
    for line in raw_file:
        if line[0] not in ['^', '!', '#', 'A']:
            stringField = line.replace('\r\n', '\n').split('\t', maxsplit= 1)
            mod_file.write(stringField[1])