#!/usr/bin/env python3.6
import os
import scipy.stats as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mlxtend.feature_selection import SequentialFeatureSelector as SFS
from sklearn.linear_model import Ridge, Lasso, ElasticNet
from sklearn.svm import SVR
import seaborn as sns
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler


def df_open(path_in, direct='n'):
    f_p, f_n = os.path.split(path_in)
    n_p, n_s = os.path.splitext(f_n)
    df1 = pd.read_csv(path_in) if n_s == '.csv' else pd.read_table(path_in)
    if direct != 'n':
        df1 = df1.T
    return df1


def scipy_corr(s1, s2, corr_mode='pearson'):
    ixs = s1.index.intersection(s2.index)
    if corr_mode == 'pearson':
        return st.pearsonr(s1[ixs], s2[ixs])
    elif corr_mode == 'kendall':
        return st.kendalltau(s1[ixs], s2[ixs])
    elif corr_mode == 'spearman':
        return st.spearmanr(s1[ixs], s2[ixs])


def fs_rfe(df_in, ss_label, f_n, tp=0):
    dfx = df_in.copy()
    if len(dfx.columns) > f_n:
        es1 = Ridge(alpha=0.1)
        es2 = Lasso(alpha=0.1)
        # es3 = SVR(C=1.0, epsilon=0.2, kernel='linear', cache_size=3000)
        ch = {0: es1, 1: es2}
        select = SFS(ch.get(tp, es1), k_features=f_n, forward=True, floating=False, verbose=1,
                     scoring='neg_mean_squared_error', cv=4, n_jobs=3)
        select.fit(dfx.values, ss_label.values)
        mask = select.k_feature_idx_
        print(mask)
        x_rfe = select.transform(dfx.values)
        m_mir_list = dfx.columns[[x for x in mask]]
        return x_rfe, m_mir_list
    else:
        f_list = dfx.columns.tolist()
        return dfx.values, f_list


def clr(tr_x, tr_y, ts_x, ts_y, path_o, fs, fn):
    train_x = tr_x.copy()
    test_x = ts_x.copy()
    clf1 = Ridge(alpha=0.1, random_state=42)
    clf2 = Lasso(alpha=0.01, max_iter=10000, random_state=42)
    clf3 = SVR(C=1.0, epsilon=0.2, kernel='linear', cache_size=3000)
    pipe1 = Pipeline([['sc', StandardScaler()], ['clf', clf1]])
    pipe2 = Pipeline([['sc', StandardScaler()], ['clf', clf2]])
    pipe3 = Pipeline([['mc', MinMaxScaler()], ['clf', clf3]])
    all_clf = [pipe1, pipe2, pipe3]
    labels = ['Ridge', 'Lasso', 'SVR']
    for clf, lb in zip(all_clf, labels):
        print('[{}]'.format(lb))
        clf.fit(train_x, tr_y)
        print('training R2:', clf.score(train_x, tr_y))
        print('testing R2', clf.score(test_x, ts_y))
        y_pre = clf.predict(X=train_x)
        yp_t = pd.Series(y_pre, index=tr_y.index)
        y_pre = clf.predict(X=test_x)
        yp_s = pd.Series(y_pre, index=ts_y.index)
        print('training cor:', scipy_corr(tr_y, yp_t))
        print('testing cor:', scipy_corr(ts_y, yp_s))
        if fn != 0:
            print('coefficients:', clf.named_steps['clf'].coef_)
            print('intercept:', clf.named_steps['clf'].intercept_)
        dfa = pd.DataFrame()
        dfa['y_testing'] = ts_y
        dfa['y_prediction'] = yp_s
        plt.rcParams["figure.figsize"] = [12, 16]
        sns.jointplot(x='y_testing', y='y_prediction', data=dfa, kind='reg')
        plt.tight_layout()
        plt.savefig(os.path.join(path_o, 'f{}_c{}_{}_testing'.format(fs, lb, fn)))
        plt.close()


def mvc(df_in, path_o, fn=0):
    train_x = df_in.iloc[:20000, 1:]
    train_y = df_in.iloc[:20000, 0]
    test_x = df_in.iloc[20000:, 1:]
    test_y = df_in.iloc[20000:, 0]
    if fn != 0:
        fn_dict = {0: 'Ridge', 1: 'LASSO', 2: 'SVR'}
        for tp in [0, 1]:
            print('Feature selection: [{}]'.format(fn_dict.get(tp)))
            tr_x, feature_l = fs_rfe(train_x, train_y, fn, tp)
            ts_x = test_x[feature_l]
            print(feature_l)
            clr(tr_x, train_y, ts_x, test_y, path_o, fn_dict.get(tp), fn)
    else:
        clr(train_x, train_y, test_x, test_y, path_o, 'No', fn)


def main():
    r_d = '.'
    f_n = 'GM_b1000_0.001.csv'
    df1 = df_open(os.path.join(r_d, f_n))
    mvc(df1, r_d, 5)


if __name__ == '__main__':
    main()
