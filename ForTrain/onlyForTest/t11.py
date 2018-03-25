#!/usr/bin/env python3.6
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from collections import defaultdict
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix
from sklearn.metrics import roc_curve
from sklearn.metrics import auc
import scipy.stats as st


a = range(30, 56, 5)
for i in a:
    print(i)


df1 = pd.DataFrame(np.random.randint(low=0, high=10, size=(4, 3)), index=['r', 'i', 'c', 'e'], columns=['m', 'i', 'r'])
ss1 = pd.Series([1, 2, 3], index=['m', 'i', 'r'])
print(df1)
print(ss1)

df2 = df1 * ss1
dfz = df1.T
dfz['cor'] = ss1
print(dfz)
print('xxx')
print(df2)
print(df2.sum(1))
df3 = df2.sort_values('m')
print(df3)
df3['x'] = range(1, 1 + len(df3.index))
print(df3)
ax = df3.plot(kind='scatter', x='x', y='m', color='b', label='mir_scores', title='Scores/phenotype')
ax.set_xlabel('rice sample index')
ax.set_ylabel('exp', color='b')
ax.tick_params('y', colors='b')
#loc = plticker.MultipleLocator(base=0.2)
#ax.xaxis.set_major_locator(loc)
ax2 = ax.twinx()
df3.plot(kind='scatter', x='x', y='r', color='r', label='phenotype', ax=ax2)
ax2.set_ylabel('sin', color='r')
ax2.tick_params('y', colors='r')
ax2.legend(loc=9)
plt.xticks(np.arange(min(df3['x']), max(df3['x'])+1, 0.2))
plt.grid()
plt.show()



def open_df(in_path, direct='n'):
    fpath, fname = os.path.split(in_path)
    fbase, fext = os.path.splitext(fname)
    df = pd.read_csv(in_path, index_col=0) if fext == '.csv' else pd.read_table(in_path, index_col=0)
    print('{}: Transpose: {}'.format(fbase, direct))
    if direct == 't':
        df = df.transpose()
    return df, fbase


d = {1: [2, 3], 2:[3, 4]}
df = pd.DataFrame(d)
print(df)


def cor_dict_get(path_in, cor_t):
    p2gp_dict, p2gm_dict = defaultdict(list), defaultdict(list)
    m2c_dict = defaultdict(dict)
    with open(path_in) as in_f:
        next(in_f)
        for line in in_f:
            lf = line.rstrip().split(',')
            if float(lf[2]) >= cor_t:
                p2gp_dict[lf[0]].append(lf[1])
                m2c_dict[lf[0]][lf[1]] = lf[2]
            elif float(lf[2]) <= -cor_t:
                p2gm_dict[lf[0]].append(lf[1])
                m2c_dict[lf[0]][lf[1]] = lf[2]
        return p2gp_dict, p2gm_dict, m2c_dict


def s_top_gt(series_input, top_n, gt=False):
    ss = series_input.copy()
    if (top_n != 0) and (top_n < len(ss)/2):
        ss.sort_values(inplace=True)
        ss = ss.iloc[range(-top_n, top_n)]
    if gt:
        ss = ss.gt(ss.mean()).astype(np.short)
    return ss


def con_m(y_true, y_pre):
    cm = confusion_matrix(y_true, y_pre)
    total = sum(sum(cm))
    print(cm)
    print(total)
    acc = (cm[0, 0] + cm[1, 1]) / total
    sensitivity1 = cm[0, 0] / (cm[0, 0] + cm[0, 1])
    specificity1 = cm[1, 1] / (cm[1, 0] + cm[1, 1])
    print('Accuracy : ', acc)
    print('Sensitivity : ', sensitivity1)
    print('Specificity : ', specificity1)
    return total, acc, sensitivity1, specificity1


def scipy_ttest_ind(s1, s2, var):
    return st.ttest_ind(s1, s2, equal_var=var)


def feature_score3(dfx, ssy_n, f_string, m2c, title_n, f_number, phe):
    dfs = dfx[f_string.split(',')]
    mc_list = [float(m2c.get(x)) for x in dfs.columns]
    ssc = pd.Series(mc_list, index=dfs.columns)
    dfz = dfs.T
    dfz['cor'] = ssc
    dfc = dfs * ssc
    df1 = pd.DataFrame()
    df1['scores'] = dfc.sum(1)
    df1[phe] = ssy_n
    df1 = df1.sort_values('scores')
    df1['index'] = range(1, 1 + len(ssy_n))
    z = np.polyfit(x=df1.loc[:, 'index'], y=df1.loc[:, phe], deg=1)
    p = np.poly1d(z)
    df1['trend'] = p(df1.loc[:, 'index'])

    plt.rcParams["figure.figsize"] = [12, 16]
    fig, ax = plt.subplots(2, 1)
    df1.plot(kind='scatter', x='index', y='scores', color='b', label='mir_scores', title='{}_{}'
                  .format(title_n, f_number), ax = ax[0])
    ax[0].axvline(x=45.5, color='g', linestyle=':')
    ax[0].text(45.5, 0.5 * (df1['scores'].min() + df1['scores'].max()), 'median score', rotation=90)
    ax[0].set_ylabel('feature scores', color='b')
    ax[0].tick_params('y')
    ax[0].grid()
    df1.plot(kind='scatter', x='index', y=phe, color='r', label='phenotype', ax=ax[1])
    df1.plot.line(x='index', y='trend', ax=ax[1])
    ax[1].axvline(x=45.5, color='g', linestyle=':')
    ax[1].text(45.5, 0.5*(df1[phe].min() + df1[phe].max()), 'median score', rotation=90)
    ax[1].set_ylabel(phe, color='r')
    ax[1].tick_params('y')
    ax[1].grid()
    plt.show()
    ss1 = df1[phe][:45]
    ss2 = df1[phe][-45:]
    print(len(ss1), len(ss2))
    print(ss1.mean(), ss2.mean())
    print(scipy_ttest_ind(ss1, ss2, False))


main_dir = 'E:/StringTemp/Project_Rice/test2'
tn = 45
dfr, n_p = open_df(os.path.join(main_dir, 'wm_all_q.csv'))
dfr = dfr.drop(['type (H)', 'waxy (H)'], axis=1)  # df2, drop all na first
dfr2 = dfr.dropna()
dfm = dfr2.iloc[:, :924]
dfp = dfr2.iloc[:, 924:]
ssp = dfp['Panicle Number (I)']
ss1 = s_top_gt(ssp, tn, True)
with open(os.path.join(main_dir, 'neg_mir.csv')) as in_f:
    for line in in_f:
        lf = line.rstrip().split(',', 2)
        if int(lf[1]) == 15:
            f_list = lf[2].split(',')
            dfx = dfm.loc[ss1.index, f_list]
            print(dfx.shape)
            clf1 = LogisticRegression(penalty='l2', C=0.001, random_state=0)
            pipe1 = Pipeline([['sc', StandardScaler()], ['clf', clf1]])
            x_train, x_test, y_train, y_test = train_test_split(dfx, ss1, test_size=0.4, random_state=1)
            log_sum = pipe1.fit(x_train, y_train)
            y_pre_tr = pipe1.predict_proba(x_train)[:, 1]   # prob include 0 and 1, only check probe of 1
            y_pre_ts = pipe1.predict_proba(x_test)[:, 1]
            n_tr, acc_tr, sen_tr, spe_tr = con_m(y_train, pipe1.predict(x_train))
            n_ts, acc_ts, sen_ts, spe_ts = con_m(y_test, pipe1.predict(x_test))
            print('\t'.join(list(map(str, ['Training', n_tr, acc_tr, sen_tr, spe_tr]))))
            print('\t'.join(list(map(str, ['Test', n_ts, acc_ts, sen_ts, spe_ts]))))
            fpr_tr, tpr_tr, thresholds_tr = roc_curve(y_true=y_train, y_score=y_pre_tr)
            fpr, tpr, thresholds = roc_curve(y_true=y_test, y_score=y_pre_ts)
            roc_auc_tr = auc(x=fpr_tr, y=tpr_tr)
            roc_auc = auc(x=fpr, y=tpr)
            print(roc_auc_tr, roc_auc)

            ssp = ssp[ss1.index]
            p2gp, p2gm, m2c = cor_dict_get('E:/StringTemp/Project_Rice/test2/wmq_corT50.csv', 0.1)
            feature_score3(dfx, ssp, lf[2], m2c.get('Panicle Number (I)'), '', 15, 'Panicle_Number_I')


dfc = pd.read_csv('E:/StringTemp/GSE50013/GSE71008_1.csv', index_col=0)
dfc = dfc[dfc.index.str.startswith('hsa-')]
dfc.rename(columns=lambda x: x.replace('GSM', 'colon'), inplace=True)
dfn = pd.read_csv('E:/StringTemp/GSE50013/GSE71008_2.csv', index_col=0)
dfn = dfn[dfn.index.str.startswith('hsa-')]
dfn.rename(columns=lambda x: x.replace('GSM', 'normal'), inplace=True)
df7 = pd.concat([dfc, dfn], axis=1)
# df7.to_csv('E:/StringTemp/GSE50013/GSE71008.csv')
