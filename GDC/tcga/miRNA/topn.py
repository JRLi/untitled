#/usr/bin/env python3.6
import pandas as pd
import os
import sys
import matplotlib.pyplot as plt
from collections import defaultdict
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFE
from sklearn.svm import SVC
from sklearn.pipeline import Pipeline
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import StandardScaler
import argparse


def args_parse():
    parser = argparse.ArgumentParser(description='no')
    parser.add_argument('file', nargs='+', help="input file")
    parser.add_argument('-m', '--mms', action="store_true", help='MinMaxScaler before feature selection')
    args = parser.parse_args()
    return args



def prepare_output_dir(output_dir):
    if os.path.exists(output_dir):
        pass
    else:
        os.mkdir(output_dir)


def df_open(path_in, direct='n'):
    f_p, f_n = os.path.split(path_in)
    n_p, n_s = os.path.splitext(f_n)
    df1 = pd.read_csv(path_in, index_col=0) if n_s == '.csv' else pd.read_table(path_in, index_col=0)
    if direct != 'n':
        df1 = df1.transpose()
    return df1, n_p


def select_r(df_in, ss_label, f_n, tp=0):
    if len(df_in.columns) > f_n:
        rfc = RandomForestClassifier(n_estimators=100, random_state=1)
        lg1 = LogisticRegression(penalty='l1', C=3, random_state=0)
        lg2 = LogisticRegression(random_state=0)
        svc = SVC(kernel='linear', probability=True, random_state=0)
        #stp = 3 if tp == 0 else 1
        stp = 1
        ch = {0: rfc, 1: lg1, 2: lg2, 3: svc}
        select = RFE(ch.get(tp, rfc), n_features_to_select=f_n, step=stp)
        select.fit(df_in, ss_label)
        mask = select.get_support()
        x_rfe = select.transform(df_in)
        m_mir_list = df_in.columns[mask]
        return x_rfe, m_mir_list, len(m_mir_list)
    else:
        f_list = df_in.columns.tolist()
        return df_in.values, f_list, len(f_list)


def mvc(dfx_in, ssy_in, base_n, mms_c):
    dfx = dfx_in.copy()
    ssy = ssy_in.copy()
    dfx1, dfx2, ssy1, ssy2 = train_test_split(dfx, ssy, test_size=0.5, random_state=0)
    print(dfx1.shape, dfx2.shape, ssy1.shape, ssy2.shape)
    fn_dict = {0: 'RandomForest', 1: 'LASSO', 2: 'Ridge', 3: 'SVC'}
    r_list = []
    f_d = defaultdict(lambda: defaultdict(dict))
    mms = MinMaxScaler()
    dfx1_m = pd.DataFrame(mms.fit_transform(dfx1), index=dfx1.index, columns=dfx1.columns) if mms_c else dfx1
    print('mms:', dfx1_m.shape)

    for f_n in range(1, 51):
        print('fn: {}'.format(f_n))
        for tp in [0, 1, 2, 3]:
            print('f_select: {}'.format(fn_dict.get(tp)))
            _, fs_list, _ = select_r(dfx1_m, ssy1, f_n, tp)
            x2 = dfx2[fs_list]
            cv = min(len(ssy2[ssy2 == 1]), len(ssy2[ssy2 == 0]))
            r_list.append(fn_dict.get(tp))

            clf1 = LogisticRegression(penalty='l2', random_state=0)
            clf2 = SVC(probability=True, random_state=0)
            clf3 = SVC(kernel='linear', probability=True, random_state=0)
            pipe1 = Pipeline([['sc', StandardScaler()], ['clf', clf1]])
            pipe2 = Pipeline([['mc', MinMaxScaler()], ['clf', clf2]])
            pipe3 = Pipeline([['sc', StandardScaler()], ['clf', clf3]])
            all_clf = [pipe1, pipe2, pipe3]
            clf_labels = ['Logistic_Regression', 'SVM_rbf', 'SVM_linear']
            for clf, label in zip(all_clf, clf_labels):
                cv_roc = cross_val_score(estimator=clf, X=x2, y=ssy2, cv=cv, scoring='roc_auc')
                f_d[label][f_n][fn_dict.get(tp)] = cv_roc.mean()
    for label in f_d.keys():
        dfg = pd.DataFrame(f_d.get(label))
        dfg.to_csv('test_{}.csv'.format(base_n))
        dfg = dfg.T
        lines = ['-', ':', '-.', '--']
        colors = ['black', 'red', 'blue', 'green']
        plt.rcParams["figure.figsize"] = [12, 9]
        for ct, clr, ls in zip(r_list, colors, lines):
            plt.plot(dfg.index, dfg[ct], color=clr, linestyle=ls, label=ct)
        plt.legend(loc='lower right')
        plt.title('Feature selection performance: {}'.format(label))
        plt.grid()
        plt.xlabel('Feature numbers')
        plt.savefig('{}_{}'.format(base_n, label))
        plt.close()


def main(argv=None):
    if argv is None:
        argv = args_parse()
        for fi in argv.file:
            df1, dfb = df_open(fi, 't')
            df1 = df1.loc[:, df1.columns.str.startswith('hsa-')]
            y_list = [0 if x.startswith('norm') else 1 for x in df1.index]
            ssy = pd.Series(y_list, index=df1.index)
            mvc(df1, ssy, dfb, argv.mms)


if __name__ == '__main__':
    sys.exit(main())
