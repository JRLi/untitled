# -*- coding: UTF-8 -*-
import pymysql.cursors

'''
import sys
accountName = str(sys.argv[1])
passWd = str(sys.argv[2])
'''
tableName = 'Standard_Diviation_15'
inputData = open("D:/Project/platelet/db.txt", "r")
for line in inputData.readlines():
    stringField = line.replace("\n", "").split("\t")
    accountName = stringField[0]
    passWd = stringField[1]
inputData.close()
conn = pymysql.connect(host='jdbc:mysql://140.120.203.44', user=accountName, password=passWd, db='siza',
                       charset='utf8mb4', cursorclass=pymysql.cursors.DictCursor)

curData = conn.cursor()
sql = 'SELECT * FROM '+tableName
curData.execute(sql)
print(curData.description)

curData.close()
conn.close()
