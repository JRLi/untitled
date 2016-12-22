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
data_dict = {'state':['ohio', 'ohio', 'ohio', 'Nevada', 'Nevada'], 'year':[2000, 2001, 2002, 2001, 2002],
             'pop':[1.5, 1.7, 3.6, 2.4, 2.9]}
frame = DataFrame(data_dict)
print(frame)
# When build df, can use 'columns' and 'index' to assign the  columns and index with order.
frame2 = DataFrame(data_dict, columns = ['year', 'state', 'pop', 'debt'], index=['one', 'two', 'three', 'four', 'five'])
print(frame2)