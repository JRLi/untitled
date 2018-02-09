from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.axes as ax
import os


def open_df(in_path, dir='n'):
    fpath, fname = os.path.split(in_path)
    fbase, fext = os.path.splitext(fname)
    df = pd.read_csv(in_path, index_col=0) if fext == '.csv' else pd.read_table(in_path, index_col=0)
    if dir == 't':
        df = df.transpose()
    return df, fbase


def sort_by_value(dict_in, rev=True):
    return sorted(dict_in, key=dict_in.get, reverse=rev)


id = 'uc001aac.4'
id2 = id[: id.find('.')]
print(id2)

a = '-'
b = '0.5'
c = float(b) if a == '+' else - float(b)
print(c)

s1 = 'ss,NaN,15,33,fsd'
lf = s1.split(',')
lf[1] = '10' if lf[1] == 'NaN' else lf[1]
print(lf)

from sklearn import linear_model
from sklearn import metrics
import numpy as np

df, df_b = open_df('E:/StringTemp/Project_Rice/nm_df129_q_i.csv')
df.drop(['type (H)', 'waxy (H)'], axis=1, inplace=True)
df_t = df.iloc[:, 924:]
df_f = df.iloc[:, :924]
print(df_t.shape)
print(' '.join(df_t.columns))
ft = 'Panicle Weight (I)'
s1 = df_t[ft]
s2 = s1.copy()
s2.sort_values(inplace=True)
ss = s2.iloc[range(-30, 30)]
print(ss)
print(ss[0], ss[-1])
ap = s2.plot(title="Panicle Weight (I)")
l_t = 'top_{}: {}'.format(30, ss[0])
l_b = 'bottom_{}: {}'.format(30, ss[-1])
ap.axhline(y=ss[0], color='g', linestyle='--', label=l_t)
ap.axhline(y=ss[-1], color='r', linestyle='--', label=l_b)
ap.legend(loc='lower right')

s1 = s1.gt(s1.mean()).astype(np.short)
lr = linear_model.LogisticRegression(penalty='l2', C=1.0)
lr.fit(df_f, s1)
preds = lr.predict_proba(df_f)
preds_1 = preds[:, 1]
fpr, tpr, _ = metrics.roc_curve(s1, preds_1)
roc_auc = metrics.auc(x=fpr, y=tpr)
plt.plot(fpr, tpr, label='{} (auc = {:.2f})'.format("Panicle Weight (I)", roc_auc))
plt.show()


df1, d_p = open_df('E:/StringTemp/Project_Rice/wm_df100k_q.csv')
df2, d_p2 = open_df('E:/StringTemp/Project_Rice/wm_df50k_q.csv')
print(df1.shape, df2.shape)
df1 = df1.dropna()
df2 = df2.dropna()
print(df1.shape, df2.shape)

x = 'g'
root = 'a' if x == 's' else 'b' if x == 'g' else 'c'
print(root)

s1 = pd.Series([1, 2, 3])
s2 = pd.Series([4, 5, 6])
s3 = s1.append(s2)
s4 = pd.concat([s1, s2])
print(s3)
print(s4)