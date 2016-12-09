import numpy as np
import pandas as pd
from sklearn import datasets, linear_model
import matplotlib.pyplot as plt

file_path = 'D:/Project/PBMC/GSE19301MCLSL'
df = pd.read_table(file_path, index_col=0)
print(df.columns)
print(df.shape)

Y = df.Disease
print(df.corr()["Disease"])
# print(Y)
# X = test_df[test_df.columns[:2]]
# print(X)
# logreg = linear_model.LogisticRegression(C=1e5)
# lg = logreg.fit(X, Y)
# print(lg)