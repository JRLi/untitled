#!/usr/bin/env python3.6
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFE
from sklearn.svm import SVC
from scipy import interpolate
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import cross_val_score
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import MinMaxScaler
r_d = 'D:\\Project\\circ_miRNA\\platform'
f_b = 'GSE113956_oral_squamous_cell_cancer.csv'
f_c = 'GSE113956_control.csv'


def open_df(path_in):
    p_d, f_n = os.path.split(path_in)
    f_p, f_s = os.path.splitext(f_n)
    df = pd.read_csv(path_in, index_col=0) if f_s == '.csv' else pd.read_table(path_in, index_col=0)
    return df


def select_r(df_in, ss_label, f_n, tp=0):
    if len(df_in.columns) > f_n > 0 and tp < 4:
        rfc = RandomForestClassifier(n_estimators=100, random_state=1)
        lg1 = LogisticRegression(penalty='l1', C=3, random_state=0)
        lg2 = LogisticRegression(penalty='l2', random_state=1)
        svc = SVC(kernel='linear', probability=True, random_state=0)
        ch = {0: rfc, 1: lg1, 2: lg2, 3: svc}
        select = RFE(ch.get(tp, lg2), n_features_to_select=f_n, step=1)
        select.fit(df_in, ss_label)
        mask = select.get_support()
        x_rfe = select.transform(df_in)
        m_list = df_in.columns[mask]
        return x_rfe, m_list
    else:
        f_list = df_in.columns.tolist()
        return df_in.values, f_list


def ml_clf(df_in, ssy_in, fn, cv):
    dfx = df_in.copy()
    ssy = ssy_in.copy()
    clf = SVC(kernel='linear', probability=True, random_state=1)
    fn_dict = {0: 'Random Forests', 1: 'Lasso Regression', 2: 'Ridge Regression', 3: 'Linear SVC', 4: 'No selection'}
    with open(os.path.join(r_d, 'mir_selected'), 'w') as out_f:
        plt.rcParams["figure.figsize"] = [12, 9]
        plt.rc('font', size=16)  # controls default text sizes
        plt.rc('axes', titlesize=16)  # fontsize of the axes title
        plt.rc('axes', labelsize=16)  # fontsize of the x and y labels
        plt.rc('xtick', labelsize=16)  # fontsize of the tick labels
        plt.rc('ytick', labelsize=16)  # fontsize of the tick labels
        plt.rc('legend', fontsize=16)  # legend fontsize
        plt.rc('figure', titlesize=16)  # fontsize of the figure title
        colors = ['k', 'r', 'g', 'b', 'orange']
        lss = ['-', '--', ':', '-.', '-.']
        for tp, color, ls in zip([4, 0, 1, 2, 3], colors, lss):
            x, f_list = select_r(dfx, ssy, fn, tp)
            cv_roc = cross_val_score(estimator=clf, X=x, y=ssy, cv=cv, scoring='roc_auc')
            out_f.write('{}\t{}\n'.format(fn_dict.get(tp), cv_roc.mean()))
            if tp in [0, 1, 2, 3]:
                out_f.write(','.join(f_list) + '\n')
            cvf = StratifiedKFold(n_splits=cv)
            y = ssy.values
            tprs, aucs= [], []
            mean_fpr = np.linspace(0, 1, 100)
            for train, test in cvf.split(x, y):
                md = clf.fit(x[train], y[train])
                pb = md.predict_proba(x[test])
                fpr_cv, tpr_cv, thresholds_cv = roc_curve(y[test], pb[:, 1])
                tprs.append(interpolate(mean_fpr, fpr_cv, tpr_cv))
                tprs[-1][0] = 0.0
                roc_auc_cv = auc(fpr_cv, tpr_cv)
                aucs.append(roc_auc_cv)
            mean_tpr = np.mean(tprs, axis=0)
            mean_tpr[-1] = 1.0
            mean_auc = auc(mean_fpr, mean_tpr)
            std_auc = np.std(aucs)
            plt.plot(mean_fpr, mean_tpr, color=color, linestyle=ls, label=r'%s (AUC = %0.2f $\pm$ %0.2f)' %
                                                         (fn_dict.get(tp), mean_auc, std_auc), lw=2, alpha=.8)
        plt.legend(loc='lower right')
        plt.plot([0, 1], [0, 1], linestyle='--', color='gray', linewidth=2)
        plt.xlim([-0.1, 1.1])
        plt.ylim([-0.1, 1.1])
        plt.title('Receiver operating characteristic of 5-Fold Linear SVC')
        plt.grid()
        plt.xlabel('False Positive Rate (1 - Specificity)')
        plt.ylabel('True Positive Rate (Sensitivity)')
        plt.tight_layout()
        plt.savefig(os.path.join(r_d, 'roc_auc'))
        plt.close()


def main():
    df_b = open_df(os.path.join(r_d, f_b))
    df_c = open_df(os.path.join(r_d, f_c))
    dfx = pd.concat([df_b, df_c], 1)
    dfx = dfx.T
    dfx = dfx.loc[(dfx != 0).any(1), (dfx != 0).any(0)]     # remove all zero data
    label_list = [1 if x in df_b.columns else 0 for x in dfx.index]
    ssy = pd.Series(label_list, index=dfx.index)
    mms = MinMaxScaler()
    dfx = pd.DataFrame(mms.fit_transform(dfx), index=dfx.index, columns=dfx.columns)
    ml_clf(dfx, ssy, 20, 5)


if __name__ == '__main__':
    main()