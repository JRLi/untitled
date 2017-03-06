import pandas as pd
from pandas import Series, DataFrame
import numpy as np
import csv

#  Series
list1 = [4, 7, -5, 3]
obj = Series(list1)
print(obj)
print(obj.values)   # Series have attribute 'values'
print(obj.index)    # Series have attribute 'index', if no assignment, default is 0, 1, 2, 3.......
# So, can change index
obj2 = Series(list1, ['d', 'b', 'a', 'c'])
print(obj2)
print(type(obj2))
print(obj2.index)
print(obj2.index[1])
print(type(obj2.index))
# So, can select values use index, and implement multiple query; like a dict, but can easy to implement multiple query
print(obj2['a'], obj[1], obj[3], obj2['d'])
print(obj2[['a', 'c', 'd']])
# and can change the values using the index
obj2['d'] = 6
obj2[['a', 'c']] = [-3, 0]
print(obj2)
# Easy to carry out complicate operator
print(obj2[obj2 > 1])
print(obj2 * 2)
print(np.exp(obj2))
# All values change, but index will keeping the same
# Series are like dictionary
print('b' in obj2, 'e' in obj2, all(['a' in obj2, 'c' in obj2]))    # all() need one iteration
data1 = {'Ohio': 35000, 'Texas': 71000, 'Oregon': 16000, 'Utah': 5000}
obj3 = Series(data1)
print(obj3)     # Attention, the default key of dict is the index, BUT have be sorted!!! The index will be sorted by
states_list = ['Oregon', 'California', 'Texas', 'Ohio']
obj4 = Series(data = data1, index = states_list)
print(obj4)     # Attention , if the index is assigned, the index will not been sorted.
# The index that is not in data_keys will be assigned NaN value
print(obj4.isnull())        # instant method!! need()
print(pd.notnull(obj4))     # pandas isnull and notnull methods
# Series adding
obj5 = obj3 + obj4
print(obj5)     # Very Important, if one index isn't in another, the value will be NaN.
# IMPORTANT, IN python, NaN's don't compare equal, but None's do.
print(None==None, np.nan==np.nan)   # True, False
a = None
b = None
print(a==b)     # True, None can be compare equal
print(obj5[obj5.isnull()])
# Process NaN
obj5[obj5.isnull()] = 0
print(obj5)     # !!!!Change the values
obj4.fillna(0)
print(obj4)     # !!! No change the value
obj4 = obj4.fillna(0)
print(obj4)     # Assign to obj4, and change the value!!!
# Series has an attribute 'name', and it's index also has the same attribute.
obj4.name='population'
obj4.index.name = 'state'
print(obj4)

# Data Frame
# Use DataFrame() to construct a df with dict
data_dict = {'state':['ohio', 'ohio', 'ohio', 'Nevada', 'Nevada', 'Nevada', 'NewYork'],
             'year':[2000, 2001, 2002, 2001, 2002, 2002, 2001],
             'pop':[1.5, 1.7, 3.6, 2.4, 2.9, 2.8, 3.9]}
frame = DataFrame(data_dict)
print(frame)
# When build df, can use 'columns' and 'index' to assign the  columns and index with order.
frame2 = DataFrame(data_dict, columns = ['year', 'state', 'pop', 'debt'],
                   index=['one', 'two', 'three', 'four', 'five', 'six', 'seven'])
print('for "for" test')
print(frame2)
print(len(frame2))

for col in frame2:  # important!!!!!!!!!!!! good fot column assess
    print(col)
    print('next')

print(frame2.iloc[1])
print(len(frame2.iloc[1]))
frame3 = frame2.transpose
print(frame3)
frame1 = frame.drop('pop', 1)
print(frame)
print(frame1)
cells = frame.columns.tolist()
print(cells)
cells.remove('state')
print('cells:', cells)
print(frame[cells])
print(frame.iloc[:,1])
df2 = frame[cells]
for i in range(0,len(df2.columns)):
    qv = np.percentile(df2.iloc[:,i], 10)
    print(qv)
    df2.iloc[:,i][df2.iloc[:,i]<= qv ] = qv
    print(df2.iloc[:,i])

df3 = pd.read_table('D://Project/circ_bicluster/CircRNABloodExp.csv', sep=',', index_col=0)
df3 = df3.iloc[0:8, 0:8]
listA = [1,3,4,6]
listB = [1,2,3]
np_array1 = np.array(listA) - 1
np_array2 = np.array(listB) - 1
print(np_array1)
print(df3.iloc[np_array1,np_array2])    # or df3.iloc[[1,3,4],[0,1,2]]
print(str(len(np_array1))+'x'+str(len(np_array2)))

# print(qv)
# filePath = 'E:/StringTemp/GDS/'
# df1 = pd.read_table(filePath + 'GDS3876.matrix', index_col=0)
# print(df1.shape)
# print(df1.index.name)
# print(df1.columns.name)
# print(df1.index, df1.columns, sep='\n')
# df1.index.name , df1.columns.name = 'genes', 'samples'
# print(df1.shape)
# print(df1.index.name)
# print(df1.columns.name)
# print(type(df1.index.name))
# print('AAAAA')
# print(df1.shape)
# print(df1.index.name)
# print(df1.columns.name)
# print(df1)
# df2 = df1.loc['ABAT'].iloc[:, 0:4]
# print(df2)
# gp = df2.groupby(df2.index)
# print(gp.mean())
# print(df1.shape[0], df1.shape[1])
# print(df1.index, df1.columns, sep='\n')
# print(df1.loc['ABAT'].iloc[:, 0:4])
# print(df1.loc['ABAT'].iloc[:, 0:4].mean(axis=0))
# print(df1.loc['ABAT'].iloc[:, 0:4].mean(axis=0).shape)
# print(df1.loc['ABAT'].iloc[:, 0:4].mean(axis=0).index)
# print(df1.index.names)
# print(df1.loc['ABAT'].iloc[:, 0:4].shape)
# a = df1.loc['ABAT'].iloc[:, 0:4].shape[0]
#
# print(df1.iloc[0:5, 0:5])
#
#
# # Numpy array test
# arr = np.array([[1, 2, 3], [4, 5, 6]])
# arr2 = np.array([[1., 2, 3], [4., 5, 6], [7, 8, 9]])
# print(arr.dtype, arr2.dtype)
# print(arr2 + 0.1)
# print(np.log2(arr + 0.1))   # use pseudo count and log2
# print(arr2 * np.array([[1], [0], [1]]))     # multiplied  with vector array
# print(arr, arr2, sep='\n')
# print(arr2.reshape(1,9))
# print(arr2)
# arr3 = arr2 * np.array([[1], [0], [1]])
# print(arr3.reshape(1, 9))
# df = pd.DataFrame(arr3, copy=True)
# arr4 = arr3.copy()
# df2 = df.copy()
# print(df)
# arr3[arr3 < 3] = 3
#
# print(df)
# print(df2)
# df2[df2 < 3] = 3
# print(df2)
#
# print(arr4)
# qv = np.percentile(arr4, 10)
# qv2 = np.percentile(arr2, 10)
# qv3 = np.percentile(arr2, 50, 1)
# print(qv, qv2, qv3)