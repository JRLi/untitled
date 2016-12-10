import numpy as np
import pandas as pd
from sklearn import datasets, linear_model
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split

# file_path = 'D:/Project/PBMC/GSE19301MCLSL'
file_path = 'D:/Rworkspace/GSE19301MCLSL'
df = pd.read_table(file_path, index_col=0)

target = 'Disease'
# print(df.corr()["Disease"])
cells = df.columns.tolist()
features = cells[1:7]   # can use for loop to do one cell iteration
print(features)
train_df, test_df = train_test_split(df, train_size=0.8, random_state=1)
print(test_df.shape, train_df.shape)
logreg = linear_model.LogisticRegression(C=1e5)
logreg.fit(train_df[features], train_df[target])
print(logreg.predict(test_df[features]).tolist(), test_df[target].tolist(), sep="\n")
print(logreg.score(train_df[features], train_df[target]))
print(logreg.score(test_df[features], test_df[target]))
print(test_df[features][:2])
print(logreg.predict_proba(test_df[features][:2]))
# print(Y)
# X = test_df[test_df.columns[:2]]
# print(X)
# logreg = linear_model.LogisticRegression(C=1e5)
# lg = logreg.fit(X, Y)
# print(lg)