#!/usr/bin/env python3.6
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn import metrics
from itertools import cycle
from numpy import interp
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import RFE
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import Imputer
from sklearn.model_selection import cross_val_score


def open_df(path_in, direct='n'):
    f_p, f_n = os.path.split(path_in)
    f_pre, f_suf = os.path.splitext(f_n)
    df1 = pd.read_csv(path_in, index_col=0) if f_suf == '.csv' else pd.read_table(path_in, index_col=0)
    if direct != 'n':
        df1 = df1.T
    return df1, f_p, f_pre


def plot_roc(fpr, tpr, roc_auc, test_series, target, des, mun=5):
    plt.figure(1)
    colors = cycle(['aqua', 'darkorange', 'cornflowerblue', 'green', 'red'])
    lw = 2
    for i, color in zip(test_series.iloc[-mun:].index, colors):
        plt.plot(fpr[i], tpr[i], color=color, lw=lw, label='{} ({:0.2f})'.format(i, roc_auc[i]))
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC of logistic regression: ' + target)
    plt.legend(loc="lower right")
    fn = '{}_{}'.format(target, des)
    plt.savefig(os.path.join('D:/Project/met/roc', fn))
    plt.gcf().clear()


def select_r(df_in, ss_label, f_n, eps):
    if len(df_in.columns) > f_n:
        select = RFE(RandomForestClassifier(n_estimators=eps, random_state=1), n_features_to_select=f_n)
        select.fit(df_in, ss_label)
        mask = select.get_support()
        x_rfe = select.transform(df_in)
        m_mir_list = df_in.columns[mask]
        return x_rfe, ','.join(m_mir_list)
    else:
        f_list = df_in.columns.tolist()
        return df_in.values, ','.join(f_list)


def logistic_single2(df_in, prefix, pen='l2', reg=1, te=0.40):
    df_a = df_in.copy()
    df_x = df_a.iloc[:, : -1]
    ss_y = df_a.iloc[:, -1]
    df_xi = df_x.fillna(df_x.mean())
    x_train, x_test, y_train, y_test = train_test_split(df_xi, ss_y, test_size=te, random_state=1)
    ny = ss_y.values
    cv = StratifiedKFold(n_splits=5)
    for i in df_x.columns:
        nX = df_xi[[i]].values
        lr = LogisticRegression(penalty=pen, C=reg)
        tprs, aucs = [], []
        mean_fpr = np.linspace(0, 1, 100)
        j = 0
        plt.rcParams["figure.figsize"] = [12, 9]
        for train, test in cv.split(nX, ny):
            md = lr.fit(nX[train], ny[train])
            pb = md.predict_proba(nX[test])
            fpr, tpr, thresholds = metrics.roc_curve(ny[test], pb[:, 1])
            tprs.append(interp(mean_fpr, fpr, tpr))
            tprs[-1][0] = 0.0
            roc_auc = metrics.auc(fpr, tpr)
            aucs.append(roc_auc)
            plt.plot(fpr, tpr, lw=1, alpha=0.3, label='ROC fold %d (AUC = %0.2f)' % (j, roc_auc))
            j += 1
        mean_tpr = np.mean(tprs, axis=0)
        mean_tpr[-1] = 1.0
        mean_auc = metrics.auc(mean_fpr, mean_tpr)
        std_auc = np.std(aucs)
        plt.plot(mean_fpr, mean_tpr, color='b', label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
                 lw=2, alpha=.8)
        std_tpr = np.std(tprs, axis=0)
        tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
        tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
        plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2, label=r'$\pm$ 1 std. dev.')
        cv_roc = cross_val_score(estimator=lr, X=df_xi[[i]], y=ss_y, cv=5, scoring='roc_auc')
        print('{}:{}'.format(i, cv_roc))
        lr.fit(x_train[[i]], y_train)   # x need to be a np.array, not a series or 1d array
        predict_y = lr.predict_proba(x_test[[i]])[:, 1]
        predict_yt = lr.predict_proba(x_train[[i]])[:, 1]
        fpr, tpr, _ = metrics.roc_curve(y_test, predict_y)
        roc_auc = metrics.auc(fpr, tpr)
        fpr_t, tpr_t, _ = metrics.roc_curve(y_train, predict_yt)
        roc_auc_t = metrics.auc(fpr_t, tpr_t)
        plt.plot(fpr, tpr, color='r', linestyle=':', label='{}-{} (auc = {:.2f})'.format(i, 'test', roc_auc), lw=2, alpha=.8)
        plt.plot(fpr_t, tpr_t, color='g', linestyle='-.', label='{}-{} (auc = {:.2f})'.format(i, 'train', roc_auc_t), lw=2, alpha=.8)
        plt.legend(loc='lower right')
        plt.plot([0, 1], [0, 1], linestyle='--', color='gray', linewidth=2)
        plt.xlim([-0.1, 1.1])
        plt.ylim([-0.1, 1.1])
        plt.title('Receiver operating characteristic of [{}]'.format(i))
        plt.grid()
        plt.xlabel('False Positive Rate (1 - Specificity)')
        plt.ylabel('True Positive Rate (Sensitivity)')
        plt.savefig(os.path.join('D:/Project/met/roc', 'ohca_{}_{}_{}c{}'.format(i, prefix, pen, reg)))
        plt.gcf().clear()


def logistic_select(df_x, ss_y, f_pre, f_number, eps, ts, reg_c=0.001):
    x = df_x.copy()
    if f_number == 0:
        feature_string = ','.join(x.columns.tolist())
        x = x.values
    else:
        x, feature_string = select_r(x, ss_y, f_number, eps)
    x_train, x_test, y_train, y_test = train_test_split(x, ss_y, test_size=ts, random_state=1)
    clf1 = LogisticRegression(penalty='l2', C=reg_c, random_state=0)
    pipe1 = Pipeline([['sc', StandardScaler()], ['clf', clf1]])
    y_pre = pipe1.fit(x_train, y_train).predict_proba(x_test)[:, 1]
    fpr, tpr, thresholds = metrics.roc_curve(y_true=y_test, y_score=y_pre)
    roc_auc = metrics.auc(x=fpr, y=tpr)
    plt.plot(fpr, tpr, color='orange', linestyle=':', label='{} (auc = {:.2f})'.format(feature_string, roc_auc))
    plt.legend(loc='lower right')
    plt.plot([0, 1], [0, 1], linestyle='--', color='gray', linewidth=2)
    plt.xlim([-0.1, 1.1])
    plt.ylim([-0.1, 1.1])
    plt.title('ROC')
    plt.grid()
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.savefig(os.path.join('D:/Project/met/roc', '{}_f{}'.format(f_pre, f_number)))
    plt.gcf().clear()


def cyb_pre(d_r):
    file_list = os.listdir(d_r)
    for f in file_list:
        if f.endswith('_2t.txt'):
            df1, f_p, df_pre = open_df(os.path.join(d_r, f))
            # print(type(df1.iloc[:, -1:]))   # -1: is df, -1 is series; iloc vs loc
            # logistic_single(df1.iloc[:, :-1], df1.iloc[:, -1], df_pre)
            logistic_select(df1.iloc[:, :-1], df1.iloc[:, -1], df_pre, 0, 100, 0.4, 1)


def sk_imputation(df_in):
    df = df_in.copy()
    imr = Imputer(missing_values='NaN', strategy='mean', axis=0)
    imr = imr.fit(df)
    arr = imr.transform(df)
    return arr


def ohca(path_in):
    df1, f_path, df_pre = open_df(path_in)
    df2 = df1[df1['DRN'] == 0]
    x1d = df1.iloc[:, np.r_[0: 19, 25]]     # np.r_ important
    x2d = df2.iloc[:, np.r_[0: 19, 25]]
    logistic_single2(x1d, 'a')
    logistic_single2(x2d, 'rm')



def main(argv=None):
    if argv is None:
        ohca('D:/Project/met/101_106_OHCA2.txt')


if __name__ == '__main__':
    sys.exit(main())
