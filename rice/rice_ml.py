#!/usr/bin/env python3.6
import os
import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
from collections import defaultdict
from sklearn import svm
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.pipeline import Pipeline
from sklearn.tree import DecisionTreeClassifier
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_curve
from sklearn.metrics import auc
from sklearn.feature_selection import RFE
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import f_regression
from sklearn.model_selection import train_test_split
from sklearn.ensemble import VotingClassifier
from sklearn.linear_model import LogisticRegression
use_message = '''
    Need Python3 and numpy, pandas, scipy, and sklearn; or Anaconda.
    Use 129 rices (113 rices if drop na) * 924 miRNAs and 17 phenotypes for machine learning.
    ps. 114 rices read counts > 50k (99 rices if drop na)
         88 rices read counts > 100k (76 rices if drop na)
    Example: python -u rice_ml.py -t 50 -ct 40 -c 0.1 -i mv -v 10 -s 0.4 -f 15
'''


def args_parse():
    parser = argparse.ArgumentParser(description=use_message)
    parser.add_argument('-m', '--mean', choices=['w', 'n'], default='n',
                        help='exp with (w) or without (n) mean, default is n')
    parser.add_argument('-r', '--read_t', choices=['a', 'f', 'h'], default='a',
                        help='total raw exp threshold of rice, "a" is all, "f" is 50K, "h" is 100K, default is a')
    parser.add_argument('-n', '--normalization', choices=['q', 'r'], default='r',
                        help='expression normalization, r is rpm, q is quantile, default is r')
    parser.add_argument('-t', '--top', type=c_int, default=50, help='Top/down -t rice; if set 0 equal all, must '
                                                                    'between 0 and 55, default is 50')
    parser.add_argument('-ct', '--cor_t', type=c_int, default=40, help='choice the correlation file type, between 0 and'
                                                                       '55, default is 40, indicate cor_t40.csv')
    parser.add_argument('-c', '--cor', type=c_float, default=0.1, help='Correlation threshold of miRNAs have between'
                                                                       'miRNA and phenotype, 0 ~ 1, default is 0.1')
    parser.add_argument('-i', '--imputation', action='store_true', help="if set, use imputation, else drop nan")
    parser.add_argument('-s', '--test_size', type=c_float, default=0.4, help="test size, 0 ~ 1, default is 0.4")
    parser.add_argument('-f', '--f_number', type=int, default=10, help='number of features, default is 10')
    subparsers = parser.add_subparsers(help='choice machine learning method', dest='command')
    parser_a = subparsers.add_parser('mv', help='method is Major voting')
    parser_a.add_argument('-v', '--cv', type=int, default=10, help='-v fold cross validation, default is 10')
    parser_b = subparsers.add_parser('svm', help='method is ANOVA linear SVM')
    parser_b.add_argument('-k', '--kernel', choices=['linear', 'poly', 'rbf', 'sigmoid'], default='linear',
                          help='Specifies the kernel type to be used in the algorithm, default is linear')
    parser_b.add_argument('-p', '--penalty', type=float, default=1.0, help='Penalty of the error term, default is 1.0')
    args = parser.parse_args()
    return args


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


def c_int(value):
    i_value = int(value)
    if i_value < 0 or i_value > 55:
        raise argparse.ArgumentTypeError("%s is an invalid int, must between 0 and 55" % value)
    return i_value


def c_float(value):
    f_value = float(value)
    if f_value < 0.0 or f_value > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]" % (f_value,))
    return f_value


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


def scipy_corr(s1, s2, corr_mode):
    ixs = s1.index.intersection(s2.index)
    if corr_mode == 'pearson':
        return st.pearsonr(s1[ixs], s2[ixs])
    elif corr_mode == 'spearman':
        return st.spearmanr(s1[ixs], s2[ixs])


def rice_corr(df_mir, df_tar, top=0, cor='pearson'):
    df_f = df_mir.copy()
    df_t = df_tar.copy()
    df_c = pd.DataFrame(columns=df_t.columns, index=df_f.columns)
    for ct in df_t.columns:
        series_t = df_t[ct]
        series_t.dropna(inplace=True)
        series_t = s_top_gt(series_t, top)
        c_list = []
        for cm in df_f.columns:
            series_c = df_f[cm]
            corr, pvl = scipy_corr(series_c, series_t, cor)
            c_list.append(corr)
        df_c[ct] = pd.Series(c_list).values
    df_c = df_c.fillna(0)
    return df_c


def percentile(df_input, per_th=10):
    df = df_input.copy()
    for i in range(0, len(df.columns)):
        qv = np.percentile(df.iloc[:, i], per_th)
        print(qv)
        df.iloc[:, i][df.iloc[:, i] < qv] = qv
    return df


def extract_rc(output, df_in, ii, tp_3):
    with open(output, 'w') as out_f:
        out_f.write('Phenotype,miRNA,{}\n'.format(tp_3))
        for r, c in zip(ii[0], ii[1]):
            phe = df_in.columns[c]
            mir = df_in.index[r]
            out_f.write('{},{},{}\n'.format(phe, mir, df_in.iloc[r, c]))


def locate(df_input, mode, threshold):
    if mode == 'p':
        ii = np.where(df_input.values <= threshold)
    else:
        ii = np.where((df_input.values >= threshold) | (df_input.values <= -threshold))     # numpy or = |
    return ii


def sort_by_value(dict_in, rev=True):
    return sorted(dict_in, key=dict_in.get, reverse=rev)    # return key list, the first element is with largest value


def plot_top(series_input, top_cor, top_ml, s_title, out_path):
    ss = series_input.copy()
    ss.sort_values(inplace=True)
    sc = ss.iloc[range(-top_cor, top_cor)]
    sm = ss.iloc[range(-top_ml, top_ml)]
    ap = ss.plot(title=s_title)
    lb_pre = 'cor&learning' if top_cor == top_ml else 'cor'
    cl_t = '{}_top_{}: {}'.format(lb_pre, top_cor, sc[0])
    cl_b = '{}_bottom_{}: {}'.format(lb_pre, top_cor, sc[-1])
    ap.axhline(y=sc[0], color='g', linestyle='--', label=cl_t)
    ap.axhline(y=sc[-1], color='g', linestyle='--', label=cl_b)
    if top_cor != top_ml:
        ml_t = 'learning_top_{}: {}'.format(top_ml, sm[0])
        ml_b = 'learning_bottom_{}: {}'.format(top_ml, sm[-1])
        ap.axhline(y=sm[0], color='r', linestyle='--', label=ml_t)
        ap.axhline(y=sm[-1], color='r', linestyle='--', label=ml_b)
    ap.legend(loc='lower right')
    fn = 'Sample_' + s_title.replace(' ', '_').replace('(', '').replace(')', '')
    plt.savefig(os.path.join(out_path, fn))
    plt.gcf().clear()


def s_top_gt(series_input, top_n, gt=False):
    ss = series_input.copy()
    if top_n != 0:
        ss.sort_values(inplace=True)
        ss = ss.iloc[range(-top_n, top_n)]
        if gt:
            ss = ss.gt(ss.mean()).astype(np.short)
    return ss


def select_r(df_in, ss_label, f_n, eps):
    if len(df_in.columns) > f_n:
        select = RFE(RandomForestClassifier(n_estimators=eps, random_state=1), n_features_to_select=f_n)
        select.fit(df_in, ss_label)
        mask = select.get_support()
        x_rfe = select.transform(df_in)
        m_mir_list = df_in.columns[mask]
        return x_rfe, ','.join(m_mir_list), len(m_mir_list)
    else:
        f_list = df_in.columns.tolist()
        return df_in.values, ','.join(f_list), len(f_list)


def anova_svm(df_in, ss_label, k_type, ts, f_n, svm_c, title_n, out_path):
    x = df_in.copy()
    x = x.values
    f_list = df_in.columns.tolist()
    fn = f_n if len(f_list) > f_n else 'all'
    x_train, x_test, y_train, y_test = train_test_split(x, ss_label, test_size=ts, random_state=1)

    anova_filter = SelectKBest(f_regression)
    clf = svm.SVC(kernel=k_type, probability=True)
    an_sv = Pipeline([('anova', anova_filter), ('svc', clf)])
    an_sv.set_params(anova__k=fn, svc__C=svm_c).fit(x_train, y_train)
    y_score = an_sv.decision_function(x_test)
    average_precision = average_precision_score(y_test, y_score)
    precision, recall, _ = precision_recall_curve(y_test, y_score)

    mask = an_sv.named_steps.anova.get_support()
    m_mir_list = df_in.columns[mask]
    f_list = ','.join(m_mir_list)
    tmp1 = '{}\t{}'.format(len(m_mir_list), f_list)

    plt.step(recall, precision, color='b', alpha=0.2, where='post')
    plt.fill_between(recall, precision, step='post', alpha=0.2, color='b')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.ylim([0.0, 1.05])
    plt.xlim([0.0, 1.0])
    plt.title('{} Precision-Recall curve: AP={:0.2f}'.format(title_n, average_precision))
    fn = 'PR_' + title_n.replace(' ', '_').replace('(', '').replace(')', '')
    plt.savefig(os.path.join(out_path, fn))
    plt.gcf().clear()

    y_pre = an_sv.predict_proba(x_test)[:, 1]
    fpr, tpr, thresholds = roc_curve(y_true=y_test, y_score=y_pre)
    roc_auc = auc(x=fpr, y=tpr)
    lf = title_n.split('_', 1)
    tmp2 = '{}\t{}\t{:.2f}'.format(lf[0], lf[1], roc_auc)
    plt.plot(fpr, tpr, color='b', linestyle='-', label='{} (auc = {:.2f})'.format(title_n, roc_auc))
    plt.legend(loc='lower right')
    plt.plot([0, 1], [0, 1], linestyle='--', color='gray', linewidth=2)
    plt.xlim([-0.1, 1.1])
    plt.ylim([-0.1, 1.1])
    plt.title(title_n)
    plt.grid()
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    fn = 'ROC_' + title_n.replace(' ', '_').replace('(', '').replace(')', '')
    plt.savefig(os.path.join(out_path, fn))
    plt.gcf().clear()
    return tmp1, tmp2


def mvc(df_x, ss_y, title_n, out_path, ts, f_number, eps, df_cv):
    x = df_x.copy()
    x, feature_string, f_l = select_r(x, ss_y, f_number, eps)
    tmp1 = '{}\t{}'.format(f_l, feature_string)
    x_train, x_test, y_train, y_test = train_test_split(x, ss_y, test_size=ts, random_state=1)
    colors = ['black', 'orange', 'blue', 'green']
    lines = [':', '--', '-.', '-']
    clf1 = LogisticRegression(penalty='l2', C=0.001, random_state=0)
    clf2 = DecisionTreeClassifier(max_depth=3, criterion='entropy', random_state=0)
    clf3 = KNeighborsClassifier(n_neighbors=5, p=2, metric='minkowski')
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
        y_pre = clf.fit(x_train, y_train).predict_proba(x_test)[:, 1]
        fpr, tpr, thresholds = roc_curve(y_true=y_test, y_score=y_pre)
        roc_auc = auc(x=fpr, y=tpr)
        plt.plot(fpr, tpr, color=clr, linestyle=ls, label='{} (auc = {:.2f})'.format(label, roc_auc))
        tmp2_list.append('{}\t{}\tROC AUC [{}]\t{:.2f}\t+/- {:.2f}\t{:.2f}'.
                         format(lf[0], lf[1], label, scores.mean(), scores.std(), roc_auc))
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
    plt.gcf().clear()
    return tmp1, tmp2


def rice_mv(df_fi, df_ti, top_n, c_t, path_c, out_p, t_size, f_n, c_v, des, ct):
    df_t = df_ti.copy()
    df_f = df_fi.copy()
    p2gp, p2gm = cor_dict_get(path_c, c_t)
    with open('mv_out/b_{}_ct{}_top{}_cor{}_t{}_f{}_cv{}'.format(des, ct, top_n, c_t, t_size, f_n, c_v), 'w') as of1, \
            open('mv_out/f_{}_ct{}_top{}_cor{}_t{}_f{}_cv{}'.format(des, ct, top_n, c_t, t_size, f_n, c_v), 'w') as of2:
        of2.write('Phenotype\tmir_cor_type\troc_method\t{}_fold_cv_auc\t{}_fold_cv_std\troc_test\n'.format(c_v, c_v))
        for t in df_t.columns:
            of1.write(t + '\n')
            of1.write('positive_cor_mir\t{}\nnegative_cor_mir\t{}\np_and_n_cor_mir\t{}\n'
                      .format(len(p2gp.get(t)), len(p2gm.get(t)), len(p2gp.get(t)) + len(p2gm.get(t))))
            plot_top(df_t[t], ct, top_n, t, out_p)
            ss1 = s_top_gt(df_t[t], top_n, True)
            of1.write('low_rice\t{}\t{}\nhigh_rice\t{}\t{}\n'
                      .format(len(ss1[ss1 == 0]), ','.join(ss1[ss1 == 0].index.tolist()), len(ss1[ss1 == 1]), ','.
                              join(ss1[ss1 == 1].index.tolist())))
            df_fp = df_f.loc[ss1.index, p2gp.get(t)]
            df_fm = df_f.loc[ss1.index, p2gm.get(t)]
            df_fa = pd.concat([df_fp, df_fm], 1)

            mir_p, roc_p = mvc(df_fp, ss1, '{}_positive'.format(t), out_p, t_size, f_n, 300, c_v)
            mir_n, roc_n = mvc(df_fm, ss1, '{}_negative'.format(t), out_p, t_size, f_n, 300, c_v)
            mir_a, roc_a = mvc(df_fa, ss1, '{}_p_and_n'.format(t), out_p, t_size, f_n, 300, c_v)
            of1.write('Positive_cor\t{}\nNegative_cor\t{}\nP_and_N_cor\t{}\n\n'.format(mir_p, mir_n, mir_a))
            of2.write('\n'.join([roc_p, roc_n, roc_a]) + '\n')


def svm_sys(df_fi, df_ti, top_n, c_t, path_c, out_p, t_size, f_n, svm_c, kernel_s, des, ct):
    df_t = df_ti.copy()
    df_f = df_fi.copy()
    p2gp, p2gm = cor_dict_get(path_c, c_t)
    with open('svm_out/b_{}_ct{}_top{}_cor{}_test{}_f{}'.format(des, ct, top_n, c_t, t_size, f_n), 'w') as of1, \
            open('svm_out/f_{}_ct{}_top{}_cor{}_test{}_f{}'.format(des, ct, top_n, c_t, t_size, f_n), 'w') as of2:
        of2.write('Phenotype\tmir_cor_type\tROC_AUC\n')
        for t in df_t.columns:
            of1.write('{}\n'.format(t))
            of1.write('positive_cor_mir\t{}\nnegative_cor_mir\t{}\np_and_n_cor_mir\t{}\n'
                      .format(len(p2gp.get(t)), len(p2gm.get(t)), len(p2gp.get(t)) + len(p2gm.get(t))))
            plot_top(df_t[t], ct, top_n, t, out_p)
            ss1 = s_top_gt(df_t[t], top_n, True)
            of1.write('low_rice\t{}\t{}\nhigh_rice\t{}\t{}\n'
                      .format(len(ss1[ss1 == 0]), ','.join(ss1[ss1 == 0].index.tolist()), len(ss1[ss1 == 1]), ','
                              .join(ss1[ss1 == 1].index.tolist())))
            df_fp = df_f.loc[ss1.index, p2gp.get(t)]
            df_fm = df_f.loc[ss1.index, p2gm.get(t)]
            df_fa = pd.concat([df_fp, df_fm], 1)

            up_r, roc_p = anova_svm(df_fp, ss1, kernel_s, t_size, f_n, svm_c, '{}_positive'.format(t), out_p)
            dn_r, roc_n = anova_svm(df_fm, ss1, kernel_s, t_size, f_n, svm_c, '{}_negative'.format(t), out_p)
            ud_r, roc_a = anova_svm(df_fa, ss1, kernel_s, t_size, f_n, svm_c, '{}_p_and_n'.format(t), out_p)

            of1.write('Positive_cor\t{}\nNegative_cor\t{}\nP_and_N_cor\t{}\n\n'.format(up_r, dn_r, ud_r))
            of2.write('\n'.join([roc_p, roc_n, roc_a]) + '\n')


def main(argv=None):
    if argv is None:
        argv = args_parse()
        prepare_output_dir('cor_file')
        prefix_f = 'nm' if argv.mean == 'n' else 'wm'
        root_f = 'all' if argv.read_t == 'a' else '50k' if argv.read_t == 'f' else '100k'
        suffix_f = argv.normalization
        extension_f = '_i.csv' if argv.imputation else '.csv'
        r_file = '_'.join([prefix_f, root_f, suffix_f]) + extension_f
        info = prefix_f + root_f + suffix_f + '_i' if argv.imputation else prefix_f + root_f + suffix_f + '_dp'
        print('\t[File]: {}'.format(r_file))

        df, df_b = open_df(r_file)
        df.drop(['type (H)', 'waxy (H)'], axis=1, inplace=True)

        if not argv.imputation:
            print('No imputation, drop all NA')
            df = df.dropna()

        df_t = df.iloc[:, 924:]
        df_f = df.iloc[:, :924]
        c_path = os.path.join('cor_file', 'cor_t{}.csv'.format(argv.cor_t))

        if not os.path.exists(c_path):
            df_c = rice_corr(df_f, df_t, argv.cor_t, 'pearson')
            extract_rc(c_path, df_c, locate(df_c, 'c', 0.01), 'correlation')

        if argv.command == 'mv':
            print('Using Major Voting')
            prepare_output_dir('mv_out')
            o_p = 'mv_out/roc_{}_ct{}_top{}_cor{}_test{}_f{}_cv{}_mv'.\
                format(info, argv.cor_t, argv.top, argv.cor, argv.test_size, argv.f_number, argv.cv)
            prepare_output_dir(o_p)
            rice_mv(df_f, df_t, argv.top, argv.cor, c_path, o_p, argv.test_size,
                    argv.f_number, argv.cv, info, argv.cor_t)
        elif argv.command == 'svm':
            print('Using Anova SVM')
            prepare_output_dir('svm_out')
            o_p = 'svm_out/roc_{}_ct{}_top{}_cor{}_test{}_f{}_svm'. \
                format(info, argv.cor_t, argv.top, argv.cor, argv.test_size, argv.f_number)
            prepare_output_dir(o_p)
            svm_sys(df_f, df_t, argv.top, argv.cor, c_path, o_p, argv.test_size,
                    argv.f_number, argv.penalty, argv.kernel, info, argv.cor_t)
        else:
            raise Usage('No command!!\nMust input mv or svm!')
    print('Done.')


if __name__ == '__main__':
    sys.exit(main())
