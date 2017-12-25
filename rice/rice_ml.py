#!/usr/bin/env python3.6
import os
import sys
import argparse
import pandas as pd
import numpy as np
from sklearn import svm
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.pipeline import Pipeline
from sklearn.tree import DecisionTreeClassifier
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_curve
from sklearn.metrics import auc
from sklearn.feature_selection import RFE
from sklearn.model_selection import train_test_split
from sklearn.ensemble import VotingClassifier
from sklearn.linear_model import LogisticRegression
import matplotlib.pyplot as plt
from collections import defaultdict


def args_parse():
    parser = argparse.ArgumentParser(description='rice miRNAs machine learning')
    parser.add_argument('-t', '--top', type=int, default=50, help='Top/down -t rice; if set 0 equal all, default is 50')
    parser.add_argument('-ct', '--cor_t', type=str, choices=['30', '40', '50'], default='40',
                        help='choice the correlation file type, default is "40", indicate cor_t40_c0.1.csv')
    parser.add_argument('-c', '--cor', type=float, default=0.1, help='Correlation threshold of miRNAs that have between'
                                                                     'miRNA and phenotype, default and min is 0.1')
    parser.add_argument('-i', '--imputation', action='store_true', help="if set, use imputation, else drop nan")
    subparsers = parser.add_subparsers(help='choice machine learning method', dest='command')
    parser_a = subparsers.add_parser('mv', help='method is Major voting')
    parser_a.add_argument('-v', '--cv', type=int, default=10, help='-v fold cross validation, default is 10')
    parser_a.add_argument('-s', '--test_size', type=float, default=0.4, help="test size, default is 0.4")
    parser_a.add_argument('-f', '--f_number', type=int, default=10, help='number of features, default is 10')
    parser_b = subparsers.add_parser('svm', help='method is loocv SVM')
    parser_b.add_argument('-k', '--kernel', choices=['linear', 'poly', 'rbf', 'sigmoid'],
                          help='Specifies the kernel type to be used in the algorithm, default is linear')
    args = parser.parse_args()
    return args


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


def prepare_output_dir(output_dir):
    if os.path.exists(output_dir):
        pass
    else:
        os.mkdir(output_dir)


def open_df(in_path, direct='n'):
    fpath, fname = os.path.split(in_path)
    fbase, fext = os.path.splitext(fname)
    df = pd.read_csv(in_path, index_col=0) if fext == '.csv' else pd.read_table(in_path, index_col=0)
    if direct == 't':
        df = df.transpose()
    return df, fbase


def loo_svm_1(df_in, ss_label, k_type):
    clf = svm.SVC(kernel=k_type)
    tp, fp, fn, tn = 0, 0, 0, 0
    for i in range(len(ss_label)):
        b = list(range(len(ss_label)))
        b.pop(i)
        print("Run{}:".format(i))
        test_x = df_in.iloc[i:i + 1]
        test_y = ss_label[i]
        train_x = df_in.iloc[b]
        train_y = ss_label[b]
        model = clf.fit(train_x, train_y)
        predict = model.predict(test_x)[0]
        print(test_y, predict)
        if test_y != 0 and predict != 0:
            tp += 1
        elif test_y == 0 and predict != 0:
            fp += 1
        elif test_y != 0 and predict == 0:
            fn += 1
        elif test_y == 0 and predict == 0:
            tn += 1
        p, r = 0, 0
        try:
            p = (tp / (tp + fp))
            r = (tp / (tp + fn))
        except:
            print('p or r have exception')
        print(", TP {}, FP {}, TN {}, FN {}, precision {:.3f}, recall {:.3f}".format(tp, fp, tn, fn, p, r))
    precision = tp / (tp + fp)
    recall = tp / (tp + fn)
    return "precision\t{:.3f}\trecall\t{:.3f}".format(precision, recall)


def cor_dict_get(path_in, cor_t):
    p2gp_dict, p2gm_dict = defaultdict(list), defaultdict(list)
    with open(path_in) as in_f:
        next(in_f)
        for line in in_f:
            lf = line.rstrip().split(',')
            if float(lf[2]) >= cor_t:
                p2gp_dict[lf[0]].append(lf[1])
            elif float(lf[2]) <= -cor_t:
                p2gm_dict[lf[0]].append(lf[1])
        return p2gp_dict, p2gm_dict


def sort_by_value(dict_in, rev=True):
    return sorted(dict_in, key=dict_in.get, reverse=rev)


def s_top_gt(series_input, top_n):
    ss = series_input.copy()
    if top_n != 0:
        ss.sort_values(inplace=True)
        ss = ss.iloc[range(-top_n, top_n)]
        ss = ss.gt(ss.mean()).astype(np.short)
    return ss


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


def mvc(df_x, ss_y, title_n, out_path, ts, f_number, eps, df_cv):
    x = df_x.copy()
    x, feature_string = select_r(x, ss_y, f_number, eps)
    tmp1 = '{}\t{}'.format(f_number, feature_string)
    x_train, x_test, y_train, y_test = train_test_split(x, ss_y, test_size=ts, random_state=1)
    colors = ['black', 'orange', 'blue', 'green']
    lines = [':', '--', '-.', '-']
    clf1 = LogisticRegression(penalty='l2', C=0.001, random_state=0)
    clf2 = DecisionTreeClassifier(max_depth=1, criterion='entropy', random_state=0)
    clf3 = KNeighborsClassifier(n_neighbors=1, p=2, metric='minkowski')
    pipe1 = Pipeline([['sc', StandardScaler()], ['clf', clf1]])
    pipe3 = Pipeline([['sc', StandardScaler()], ['clf', clf3]])
    mv_clf = VotingClassifier([('lr', pipe1), ('dt', clf2), ('knn', pipe3)], voting='soft')
    all_clf = [pipe1, clf2, pipe3, mv_clf]
    clf_labels = ['Logistic Regression', 'Decision Tree', 'KNN', 'Majority Voting']
    tmp2_list = []
    cv = x_train.shape[0] if df_cv > x_train.shape[0] else df_cv
    for clf, label, clr, ls in zip(all_clf, clf_labels, colors, lines):
        scores = cross_val_score(estimator=clf, X=x_train, y=y_train, cv=cv, scoring='roc_auc')
        lf = title_n.split('_', 1)
        tmp2_list.append('{}\t{}\tROC AUC [{}]\t{:.2f}\t+/- {:.2f}'.format(lf[0], lf[1], label, scores.mean(), scores.std()))
        y_pre = clf.fit(x_train, y_train).predict_proba(x_test)[:, 1]
        fpr, tpr, thresholds = roc_curve(y_true=y_test, y_score=y_pre)
        roc_auc = auc(x=fpr, y=tpr)
        plt.plot(fpr, tpr, color=clr, linestyle=ls, label='{} (auc = {:.2f})'.format(label, roc_auc))
    tmp2 = '\n'.join(tmp2_list)
    plt.legend(loc='lower right')
    plt.plot([0, 1], [0, 1], linestyle='--', color='gray', linewidth=2)
    plt.xlim([-0.1, 1.1])
    plt.ylim([-0.1, 1.1])
    plt.title(title_n)
    plt.grid()
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    fn = title_n.replace(' ', '_').replace('(', '').replace(')', '')
    plt.savefig(os.path.join(out_path, fn))
    #plt.show()
    plt.gcf().clear()
    return tmp1, tmp2


def rice_mv(df_fi, df_ti, top_n, c_t, path_c, out_p, t_size, f_n, c_v, des, ct):
    df_t = df_ti.copy()
    df_f = df_fi.copy()
    p2gp, p2gm = cor_dict_get(path_c, c_t)
    with open('b_{}_ct{}_top{}_cor{}_test{}_f{}_cv{}'.format(des, ct, top_n, c_t, t_size, f_n, c_v), 'w') as out_f, \
            open('f_{}_ct{}_top{}_cor{}_test{}_f{}_cv{}'.format(des, ct, top_n, c_t, t_size, f_n, c_v), 'w') as out_f2:
        out_f2.write('Phenotype\tmir_cor_type\troc_method\t{}_fold_cv_auc\t{}_fold_cv_std\n'.format(c_v, c_v))
        for t in df_t.columns:
            out_f.write(t + '\n')
            out_f.write('positive_cor_mir\t{}\nnegative_cor_mir\t{}\np_and_n_cor_mir\t{}\n'.
                        format(len(p2gp.get(t)), len(p2gm.get(t)), len(p2gp.get(t)) + len(p2gm.get(t))))
            ss1 = s_top_gt(df_t[t], top_n)
            out_f.write('low_rice\t{}\t{}\nhigh_rice\t{}\t{}\n'.
                        format(len(ss1[ss1 == 0]), ','.join(ss1[ss1 == 0].index.tolist()), len(ss1[ss1 == 1]),
                               ','.join(ss1[ss1 == 1].index.tolist())))
            df_fp = df_f.loc[ss1.index, p2gp.get(t)]
            df_fm = df_f.loc[ss1.index, p2gm.get(t)]
            df_fa = pd.concat([df_fp, df_fm], 1)
            mir_p, roc_p = mvc(df_fp, ss1, '{}_positive'.format(t), out_p, t_size, f_n, 300, c_v)
            mir_n, roc_n = mvc(df_fm, ss1, '{}_negative'.format(t), out_p, t_size, f_n, 300, c_v)
            mir_a, roc_a = mvc(df_fa, ss1, '{}_p_and_n'.format(t), out_p, t_size, f_n, 300, c_v)
            out_f.write('Positive_cor\t{}\nNegative_cor\t{}\nP_and_N_cor\t{}\n\n'.format(mir_p, mir_n, mir_a))
            out_f2.write('\n'.join([roc_p, roc_n, roc_a]) + '\n')


def main(argv=None):
    if argv is None:
        argv = args_parse()
        c_file = 'cor_t{}_c0.1.csv'.format(argv.cor_t)
        r_file = 'nm_df129_imputation.csv' if argv.imputation else 'nm_df129.csv'
        df, df_b = open_df(r_file)
        df.drop(['type (H)', 'waxy (H)'], axis=1, inplace=True)
        if not argv.imputation:
            print('No imputation, drop all NA')
            df = df.dropna()
        i = 'imputation' if argv.imputation else 'drop_na'
        df_t = df.iloc[:, 924:]
        df_f = df.iloc[:, :924]
        if argv.command == 'mv':
            print('Using Major Voting')
            o_p = 'roc_{}_ct{}_top{}_cor{}_test{}_f{}_cv{}_p'.\
                format(i, argv.cor_t, argv.top, argv.cor, argv.test_size, argv.f_number, argv.cv)
            prepare_output_dir(o_p)
            rice_mv(df_f, df_t, argv.top, argv.cor, c_file, o_p, argv.test_size, argv.f_number, argv.cv, i, argv.cor_t)
        elif argv.command == 'svm':
            print('Using leave one out cross validation SVM')
        else:
            raise Usage('No command!!\nMust input mv or svm!')
    print('Done.')


if __name__ == '__main__':
    sys.exit(main())