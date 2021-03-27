from sklearn import svm
from sklearn.model_selection import train_test_split
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import f_regression
from sklearn.pipeline import Pipeline
import pandas as pd
import os
import numpy as np


def open_df(in_path, direct='n'):
    fpath, fname = os.path.split(in_path)
    fbase, fext = os.path.splitext(fname)
    df = pd.read_csv(in_path, index_col=0) if fext == '.csv' else pd.read_table(in_path, index_col=0)
    if direct == 't':
        df = df.transpose()
    return df, fbase


'''
df, dfb = open_df('E:/StringTemp/Project_Rice/wm_df129_q_i.csv')
df_t = df.iloc[:, 924:]
df_f = df.iloc[:, :20]
print(df_f.shape)
x_train, x_test, y_train, y_test = train_test_split(df_f.values, df_t['Panicle Length (H)'], test_size=0.4, random_state=1)
anova_filter = SelectKBest(f_regression, k=5)
clf = svm.SVC(kernel='linear')
anova_svm = Pipeline([('anova', anova_filter), ('svc', clf)])
anova_svm.set_params(anova__k='all', svc__C=.1).fit(x_train, y_train)
prediction = anova_svm.predict(x_train)
print(anova_svm.score(x_train, y_train))
'''
df = pd.read_csv('C:/Users/fanic/Downloads/tmp/test.csv', index_col=0)
df2 = pd.read_csv('C:/Users/fanic/Downloads/tmp/test2.csv', index_col=0)
df3 = pd.read_csv('C:/Users/fanic/Downloads/tmp/test3.csv', index_col=0)
print(df)
print(df2)
idx = df.index.intersection(df2.index)
idx2 = df2.index.intersection(df.index)
print(idx)
print(idx2)
df4 = df2.loc[idx]
print(df4)
idx = df3.index.intersection(df2.index)
df5 = df2.loc[idx]
print(df5)
df5.index = df5.y
print(df5)
print((df4.x <= df.A) & (df.A < df4.y))
#df3 = pd.DataFrame()
#df3['s1'] = (df.D <= df.A) & (df.A < df.E)
#df3['s2'] = (df.D <= df.B) & (df.B < df.E)
#df3['s3'] = (df.D <= df['A.1']) & (df['A.1'] < df.E)

#print(df2.sum())


