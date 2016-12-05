# -*- coding: UTF-8 -*-
import pymysql.cursors
import csv
'''
import sys
accountName = str(sys.argv[1])
passWd = str(sys.argv[2])
'''
tableName = 'Standard_Diviation_15'
dbName = 'siza'

with open("D:/Project/platelet/db.txt", "r") as inputData:
    for line in inputData:
        stringField = line.replace("\n", "").split("\t")
        accountName = stringField[0]
        passWd = stringField[1]

conn = pymysql.connect(host='140.120.203.44', user=accountName, password=passWd, db=dbName,
                       charset='utf8', cursorclass=pymysql.cursors.Cursor)  # if dict out: pymysql.cursors.DictCursor
sql = 'SELECT * FROM ' + tableName

curData = conn.cursor()
curData.execute(sql)

print(curData.description)
field_names = [i[0] for i in curData.description]   # if need to get column name
print(field_names)
print(type(curData))    # tuple

result = curData.fetchall()
print('result type:', type(result))
print(result[1:2:])
testA = result[1]
print(testA, type(testA))
print(testA[1])
with open('E:/StringTemp/sqlTest2.csv', 'w') as fp:
    for i in field_names:       # if need to write column name
        fp.write(i + ",")
    fp.write("\n")
    myFile = csv.writer(fp, lineterminator='\n')
    myFile.writerows(result)

curData.close()
conn.close()

