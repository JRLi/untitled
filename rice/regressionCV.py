#!/usr/bin/env python3.6
import os
import pandas as pd
import numpy as np
from sklearn import svm
from sklearn.datasets import load_breast_cancer
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.pipeline import Pipeline
from sklearn.tree import DecisionTreeClassifier
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.metrics import roc_curve
from sklearn.metrics import auc
from sklearn.feature_selection import SelectPercentile, SelectFromModel, RFE
from sklearn.model_selection import train_test_split
from sklearn.ensemble import VotingClassifier
from sklearn.linear_model import LogisticRegression
import matplotlib.pyplot as plt
from collections import defaultdict


def open_df(in_path, direct='n'):
    fpath, fname = os.path.split(in_path)
    fbase, fext = os.path.splitext(fname)
    df = pd.read_csv(in_path, index_col=0) if fext == '.csv' else pd.read_table(in_path, index_col=0)
    if direct == 't':
        df = df.transpose()
    return df, fbase, fpath


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


def s_top_gt(series_input, top_n):
    ss = series_input.copy()
    if top_n != 0:
        ss.sort_values(inplace=True)
        ss = ss.iloc[range(-top_n, top_n)]
        ss = ss.gt(ss.mean()).astype(np.short)
    return ss


def loo_svm_2(df_in, ss_label):
    runs = df_in.copy()
    clf = svm.SVC(kernel='linear')
    tp, fp, fn, tn = 0, 0, 0, 0
    runs.index = ss_label
    runs.columns = range(len(runs.columns))
    for i in range(len(runs)):
        print("Run {}".format(i))
        run = runs[0:1]
        label = run.index[0]
        runs = runs[1:]  # change runs
        labels = runs.index
        train = runs.copy()
        test = run.copy()
        model = clf.fit(np.array(train), np.array(labels))
        predict = model.predict(np.array(test))[0]
        print(label, predict, sep='\t')
        runs = runs.append(run)
        if label != 0 and predict != 0:
            tp += 1
        elif label == 0 and predict != 0:
            fp += 1
        elif label != 0 and predict == 0:
            fn += 1
        elif label == 0 and predict == 0:
            tn += 1
        p = 0
        r = 0
        try:
            p = (tp / (tp + fp))
        except:
            pass
        try:
            r = (tp / (tp + fn))
        except:
            pass
        print(", TP {}, FP {}, TN {}, FN {}, precision {:.3f}, recall {:.3f}".format(tp, fp, tn, fn, p, r))
    precision = tp / (tp + fp)
    recall = tp / (tp + fn)
    return "precision\t{:.3f}\trecall\t{:.3f}".format(precision, recall)


def loo_svm_1(df_in, ss_label):
    clf = svm.SVC(kernel='linear')
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
        p = 0
        r = 0
        try:
            p = (tp / (tp + fp))
        except:
            pass
        try:
            r = (tp / (tp + fn))
        except:
            pass
        print(", TP {}, FP {}, TN {}, FN {}, precision {:.3f}, recall {:.3f}".format(tp, fp, tn, fn, p, r))
    precision = tp / (tp + fp)
    recall = tp / (tp + fn)
    return "precision\t{:.3f}\trecall\t{:.3f}".format(precision, recall)


def loo_svm_sys(df_fi, df_ti, top_n, c_t, path_c, out_p):
    with open(os.path.join(out_p, 'loosv_svm_result'), 'w') as out_f:
        out_f.write('features_shape:\t{}\ntarget_shape:\t{}\n'.format(df_fi.shape, df_ti.shape))
        out_f.write('[loosv svm]\ncor_threshold:\t{}\ntop_target:\t{}\n\n'.format(c_t, top_n))
        df_t = df_ti.copy()
        df_f = df_fi.copy()
        p2gp, p2gm = cor_dict_get(path_c, c_t)
        for t in df_t.columns:
            #if t.startswith('Panicle Number (I)'):
            out_f.write('[{}]\n'.format(t))
            ss1 = s_top_gt(df_t[t], top_n)
            out_f.write('low\t{}\t{}\nhigh\t{}\t{}\n'.format(len(ss1[ss1 == 0]), ','.join(ss1[ss1 == 0].index.tolist()),
                                                             len(ss1[ss1 == 1]), ','.join(ss1[ss1 == 1].index.tolist())))
            df_fp = df_f.loc[ss1.index, p2gp.get(t)]
            df_fm = df_f.loc[ss1.index, p2gm.get(t)]
            print(df_fp.shape, df_fm.shape)
            df_fa = pd.concat([df_fp, df_fm], 1)
            up_r = loo_svm_2(df_fp, ss1)
            dn_r = loo_svm_2(df_fm, ss1)
            ud_r = loo_svm_2(df_fa, ss1)
            out_f.write('cor >= {}\t{}\t{}\n{}\n'.format(c_t, len(p2gp.get(t)), ','.join(p2gp.get(t)), up_r))
            out_f.write('cor <= -{}\t{}\t{}\n{}\n'.format(c_t, len(p2gm.get(t)), ','.join(p2gm.get(t)), dn_r))
            out_f.write('cor+-{}\t{}\n{}\n'.format(c_t, len(p2gp.get(t)) + len(p2gm.get(t)), ud_r))


def rice_mv(df_fi, df_ti, top_n, c_t, path_c, out_p, t_size=0.5):
    df_t = df_ti.copy()
    df_f = df_fi.copy()
    p2gp, p2gm = cor_dict_get(path_c, c_t)
    print(len(p2gp), len(p2gm))
    for k, v in p2gp.items():
        print(k)
        print('p: {}\tm: {}\tp_and_m: {}'.format(len(v), len(p2gm.get(k)), len(v) + len(p2gm.get(k))))
    for t in df_t.columns:
        ss1 = s_top_gt(df_t[t], top_n)
        print(t)
        print('low \t{}\t{}\nhigh\t{}\t{}'.format(len(ss1[ss1 == 0]), ','.join(ss1[ss1 == 0].index.tolist()),
                                                         len(ss1[ss1 == 1]), ','.join(ss1[ss1 == 1].index.tolist())))
        df_fp = df_f.loc[ss1.index, p2gp.get(t)]
        df_fm = df_f.loc[ss1.index, p2gm.get(t)]
        df_fa = pd.concat([df_fp, df_fm], 1)
        print(df_fp.shape, df_fm.shape, df_fa.shape)
        #np_p, m_p = select_r(df_fp, ss1, 20)
        np_m, m_m = select_r(df_fm, ss1, 20)
        #np_a, m_a = select_r(df_fa, ss1, 20)
        #p_mir_list = df_fp.columns[m_p]
        m_mir_list = df_fm.columns[m_m]
        #a_mir_list = df_fa.columns[m_a]
        #mvc(np_p, ss1, '{}_positive'.format(t))
        mvc(np_m, ss1, '{}_negative'.format(t))
        #mvc(np_a, ss1, '{}_p_and_n'.format(t))


def select_r(df_in, ss_label, f_n=15, eps=100):
    if len(df_in.columns) > f_n:
        select = RFE(RandomForestClassifier(n_estimators=eps, random_state=1), n_features_to_select=f_n)
        select.fit(df_in, ss_label)
        mask = select.get_support()
        x_rfe = select.transform(df_in)
        return x_rfe, mask
    else:
        return df_in, None


def mvc(x, y, title_n, df_cv=10):
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
    print('10-fold cross validation:\n')
    for clf, label, clr, ls in zip(all_clf, clf_labels, colors, lines):
        scores = cross_val_score(estimator=clf, X=x, y=y, cv=df_cv, scoring='roc_auc')
        print('ROC AUC: {:.2f} (+/- {:.2f}) [{}]'.format(scores.mean(), scores.std(), label))
        #y_pre = clf.fit(x, y).predict_proba(x)[:, 1]
        y_pre = clf.predict_proba(x)[:, 1]
        fpr, tpr, thresholds = roc_curve(y_true=y, y_score=y_pre)
        roc_auc = auc(x=fpr, y=tpr)
        plt.plot(fpr, tpr, color=clr, linestyle=ls, label='{} (auc = {:.2f})'.format(label, roc_auc))
    plt.legend(loc='lower right')
    plt.plot([0, 1], [0, 1], linestyle='--', color='gray', linewidth=2)
    plt.xlim([-0.1, 1.1])
    plt.ylim([-0.1, 1.1])
    plt.title(title_n)
    plt.grid()
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.show()


def rice():
    #r_p = '.'
    r_p = 'E:/StringTemp/Project_Rice/'
    f_i = 'nm_df129_imputation.csv'
    df, df_b, df_p = open_df(os.path.join(r_p, f_i))
    df.drop(['type (H)', 'waxy (H)'], axis=1, inplace=True)
    df_t = df.iloc[:, 924:]
    df_f = df.iloc[:, :924]
    #loo_svm_sys(df_f, df_t, 30, 0.20, 'cor_t30_c0.1.csv', df_p)
    rice_mv(df_f, df_t, 50, 0.10, os.path.join(r_p, 'cor_t30_c0.1.csv'), df_p)


def sel_per(X_train, y_train, X_test, y_test):
    sel_per = SelectPercentile(percentile=50)   # use f_classif and 50% percentile
    sel_per.fit(X_train, y_train)
    X_train_selected = sel_per.transform(X_train)  # select X train
    print('X_train_shape: {}'.format(X_train.shape))
    print('X_train_selected.shape: {}'.format(X_train_selected.shape))
    mask = sel_per.get_support()
    X_test_selected = sel_per.transform(X_test)
    lr = LogisticRegression()
    lr.fit(X_train, y_train)
    print('LR score with all features: {:.3f}'.format(lr.score(X_test, y_test)))
    lr.fit(X_train_selected, y_train)
    print('LR score with selected features: {:.3f}'.format(lr.score(X_test_selected, y_test)))
    return mask


def sel_mod(X_train, y_train, X_test, y_test):
    select = SelectFromModel(RandomForestClassifier(n_estimators=100, random_state=42), threshold='median')
    select.fit(X_train, y_train)
    X_train_l1 = select.transform(X_train)
    print('X_train_shape: {}'.format(X_train.shape))
    print('X_train_l1.shape: {}'.format(X_train_l1.shape))
    mask = select.get_support()
    X_test_l1 = select.transform(X_test)
    score = LogisticRegression().fit(X_train_l1, y_train).score(X_test_l1, y_test)
    print('LR test score with mod selected features: {:.3f}'.format(score))
    return mask


def sel_rfe(X_train, y_train, X_test, y_test):
    select = RFE(RandomForestClassifier(n_estimators=100, random_state=42), n_features_to_select=40)
    select.fit(X_train, y_train)
    mask = select.get_support()
    X_train_rfe = select.transform(X_train)
    X_test_rfe = select.transform(X_test)
    score = LogisticRegression().fit(X_train_rfe, y_train).score(X_test_rfe, y_test)
    print('LR test score with RFE selected features: {:.3f}'.format(score))
    return mask


def breast():
    cancer = load_breast_cancer()
    print(cancer.target)
    print(cancer.target_names)
    rng = np.random.RandomState(42)
    # add noises
    noise = rng.normal(size=(len(cancer.data), 50))
    X_w_noise = np.hstack([cancer.data, noise])
    X_train, X_test, y_train, y_test = train_test_split(X_w_noise, cancer.target, random_state=0, test_size=.5)
    #mask = sel_per(X_train, y_train, X_test, y_test)
    #mask = sel_mod(X_train, y_train, X_test, y_test)
    mask = sel_rfe(X_train, y_train, X_test, y_test)
    print(mask)
    plt.matshow(mask.reshape(1, -1), cmap='gray_r')
    plt.xlabel('sample_index')
    plt.show()


def main():
    rice()
    #breast()


if __name__ == '__main__':
    main()
