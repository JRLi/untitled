import numpy as np
import pandas as pd
from pandas import Series, DataFrame
from scipy import stats
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt

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


file_path = 'D:/Project/PBMC/finalCLS/H1_GSE22886/GSE3365_GPL96.txt_p1000_H1_final'
input = pd.read_table(file_path, index_col=0, nrows = 10)
print(input.columns.values)
print(list(input.index))
print(type(input))
print(input)


file_path = 'D:/Project/python_learning/games.csv'
# games = pd.read_csv(file_path, index_col= 0)
games = pd.read_csv(file_path)  # no use index_col=0 because the first column isn't row name
print(games.columns)    # Print the names of the columns in games.
print(games.shape)

# Make a histogram of all the ratings in the average_rating column.
# plt.hist(games.average_rating)
# plt.hist(games['average_rating'])   # same with up described
# plt.show()    # Show the plot.

# Print the first row of all the games with zero scores.
# The .iloc method on dataframes allows us to index by position.
print(games[games["average_rating"] == 0].iloc[0])
#  print(games[games["average_rating"] > 0].iloc[0])

# Remove any rows without user reviews.
games = games[games["users_rated"] > 0]
# Remove any rows with missing values.
games = games.dropna(axis=0)

# Import the kmeans clustering model.
# Initialize the model with 2 parameters -- number of clusters and random state.
kmeans_model = KMeans(n_clusters=5, random_state=1)
# Get only the numeric columns from games.
good_columns = games._get_numeric_data()
# Fit the model using the good columns.
kmeans_model.fit(good_columns)
# Get the cluster assignment.
labels = kmeans_model.labels_

