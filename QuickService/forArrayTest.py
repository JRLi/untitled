import numpy as np
import pandas as pd
from pandas import Series, DataFrame
from scipy import stats

data1 = [6, 7.5, 8, -1, 0, 2]
arr1 = np.array(data1)
print(data1)
print(arr1)

arang1 = np.arange(1, 3, 0.5)
print(arang1)
for i in arang1:
    print(i)
# eye matrix
eye0 = np.eye(2, 4, 0)
eye1 = np.eye(2, 4, 1)
eyeS4 = np.eye(4, 4)
print(eye0)
print(eye1)
print(eyeS4)

prng = np.random.RandomState()      # random seed
rng1 = prng.permutation(10)
prng = np.random.RandomState()      # reset random seed will NOT give same number
rng2 = prng.permutation(10)
prng = np.random.RandomState(0)     # set seed 0
rng3 = prng.permutation(10)
prng = np.random.RandomState(0)     # reset seed 0, give same number
rng4 = prng.permutation(10)
rng5 = prng.permutation(10).reshape((2, 5))     # reshape to 2 by 3 matrix
print(type(prng), type(rng1))
print(rng1, rng2, rng3, rng4, sep="\n")
print(type(rng5), rng5, sep="\n")

rvs1 = stats.norm.rvs(loc=10, scale=20, size=100)
rvs2 = stats.norm.rvs(loc=0, scale=1, size=80)
tt_ind_result = stats.ttest_ind(rvs1, rvs2)
print(tt_ind_result)
print(tt_ind_result.statistic)
print(tt_ind_result.pvalue)
print()
# pd test: Series
obj = Series([4, 7, -3, 5])
print(obj)
print(obj.values)
print(obj.index)
obj2 = Series([4, 7, -5, 3], index=['a', 'b', 'c', 'd'])
print(obj2)
print(obj2.index)
print(">0 test:", obj2[obj2 > 0])
print("log2 test:", np.square(obj2))
print("key in test:", 'b' in obj2)
obj2.index = ['z', 'y', 'x', 'w']   # change index
print("index assign test:", obj2)


