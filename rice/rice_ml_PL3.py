#!/usr/bin/env python3.6
import pandas as pd
import numpy as np
import scipy.stats as st
import os
import math
import matplotlib.pyplot as plt
from numpy import interp
from collections import defaultdict
from sklearn import metrics
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFE
from mlxtend.feature_selection import SequentialFeatureSelector as SFS
from sklearn.svm import SVC
from sklearn.pipeline import Pipeline
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import roc_curve
from sklearn.metrics import auc
from sklearn.metrics import confusion_matrix
rs_stat = 25


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


def rice_corr(df_mir, df_tar, top=0, cor='pearson'):
    df_f = df_mir.copy()
    df_t = df_tar.copy()
    df_c = pd.DataFrame(columns=df_t.columns, index=df_f.columns)
    df_p = pd.DataFrame(columns=df_t.columns, index=df_f.columns)
    for ct in df_t.columns:
        series_t = df_t[ct]
        series_t.dropna(inplace=True)
        series_t = s_top_gt(series_t, top)
        c_list, p_list = [], []
        for cm in df_f.columns:
            series_c = df_f[cm]
            corr, pvl = scipy_corr(series_c, series_t, cor)
            c_list.append(corr)
            p_list.append(pvl)
        df_c[ct] = pd.Series(c_list).values
        df_p[ct] = np.array(p_list)
    df_c = df_c.fillna(0)
    df_p = df_p.fillna(1)
    return df_c, df_p


def extract_rc(output, df_in, ii, tp):
    with open(output, 'w') as out_f:
        out_f.write('Phenotype,miRNA,{}\n'.format(tp))
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


def scipy_corr(s1, s2, corr_mode):
    ixs = s1.index.intersection(s2.index)
    if corr_mode == 'pearson':
        return st.pearsonr(s1[ixs], s2[ixs])
    elif corr_mode == 'spearman':
        return st.spearmanr(s1[ixs], s2[ixs])


def s_top_gt(series_input, top_n, gt=False):
    ss = series_input.copy()
    if (top_n != 0) and (top_n < len(ss)/2):
        ss.sort_values(inplace=True)
        ss = ss.iloc[range(-top_n, top_n)]
    if gt:
        ss = ss.gt(ss.mean()).astype(np.short)
        # ss = ss.gt(11).astype(np.short)
    return ss


def con_m(y_true, y_pre):
    cm = confusion_matrix(y_true, y_pre)
    total = sum(sum(cm))
    print(cm)
    print(total)
    acc = (cm[0, 0] + cm[1, 1]) / total
    sensitivity1 = cm[1, 1] / (cm[1, 1] + cm[1, 0])
    specificity1 = cm[0, 0] / (cm[0, 0] + cm[0, 1])
    precision1 = cm[1, 1] / (cm[1, 1] + cm[0, 1])
    print('Accuracy : ', acc)
    print('Sensitivity (recall) : ', sensitivity1)
    print('Specificity : ', specificity1)
    return total, acc, sensitivity1, specificity1, precision1


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


def p_dict_get(path_in, p_t):
    m2p_dict = defaultdict(dict)
    with open(path_in) as in_f:
        next(in_f)
        for line in in_f:
            lf = line.rstrip().split(',')
            if float(lf[2]) <= p_t:
                m2p_dict[lf[0]][lf[1]] = lf[2]
        return m2p_dict


def select_r2(df_in, ss_label, f_n, eps):
    dfx = df_in.copy()
    if len(dfx.columns) > f_n:
        select = SFS(RandomForestClassifier(n_estimators=eps, random_state=1), k_features=f_n, forward=True,
                     floating=False,scoring='accuracy', cv=4, n_jobs=3)
        select.fit(dfx.values, ss_label.values)
        mask = select.k_feature_idx_
        x_sfs = select.transform(dfx.values)
        m_mir_list = dfx.columns[[x for x in mask]]
        return x_sfs, ','.join(m_mir_list), len(m_mir_list)
    else:
        f_list = dfx.columns.tolist()
        return dfx.values, ','.join(f_list), len(f_list)


def select_r(df_in, ss_label, f_n, eps):
    dfx = df_in.copy()
    if len(dfx.columns) > f_n:
        select = RFE(RandomForestClassifier(n_estimators=eps, random_state=1), n_features_to_select=f_n)
        select.fit(dfx, ss_label)
        mask = select.get_support()
        x_rfe = select.transform(dfx)
        m_mir_list = dfx.columns[mask]
        return x_rfe, ','.join(m_mir_list), len(m_mir_list)
    else:
        f_list = dfx.columns.tolist()
        return dfx.values, ','.join(f_list), len(f_list)


def fn_check_plot(df_in, path_o, title_n):
    l_path = os.path.join(path_o, 'roc_lg')
    ml_path = os.path.join(path_o, 'roc_svmL')
    prepare_output_dir(l_path)
    prepare_output_dir(ml_path)
    df1 = df_in.copy()
    for t, o_dir in zip(['Logistic Regression', 'SVM Linear'], [l_path, ml_path]):
        lines = [':', '-.', '--']
        colors = ['red', 'blue', 'green']
        col_l = ['lr_cv10_roc', 'lr_train_roc', 'lr_test_roc'] if t == 'Logistic Regression' \
            else ['svl_cv10_roc', 'svl_train_roc', 'svl_test_roc']
        for ct, ls, clr in zip(col_l, lines, colors):
            plt.plot(df1['mir_number'], df1[ct], color=clr, linestyle=ls, label='{}_{}'.format(t, ct.split('_', 1)[1]))
        plt.legend(loc='lower right')
        plt.title(title_n)
        plt.grid()
        plt.xticks(np.arange(0, max(df1['mir_number']) + 1, 5))
        plt.xlabel('Feature numbers')
        f_name = title_n.replace(' ', '_').replace('(', '').replace(')', '')
        plt.tight_layout()
        plt.savefig(os.path.join(o_dir, '{}_{}'.format(f_name, t.replace(' ', '_'))))
        plt.close()


def top_n_test(df_xi, ss_y, title_n, out_path, ts, f_number, eps, df_cv):
    plt.rcParams["figure.figsize"] = [22, 16]
    fig, ax = plt.subplots(3, 3)
    i, j = 0, 0
    list_45l, list_50l = [], []
    for tn in range(30, 56, 5):
        print('[top {}]'.format(tn))
        ss1 = s_top_gt(ss_y, tn, True)
        dfx = df_xi.loc[ss1.index, :]
        fn_list, l_tr_list, l_cv_list, l_ts_list = [], [], [], []
        ml_tr_list, ml_cv_list, ml_ts_list = [], [], []
        print('Now running F_number:')
        for fn in range(1, f_number + 1):
            print(fn, end=',')
            fn_list.append(fn)
            x, feature_string, f_l = select_r(dfx, ss1, fn, eps)
            x_train, x_test, y_train, y_test = train_test_split(x, ss1, test_size=ts, random_state=rs_stat)
            clf1 = LogisticRegression(penalty='l2', C=0.001, random_state=0)
            clf3 = SVC(kernel='linear', probability=True, random_state=0)
            pipe1 = Pipeline([['sc', StandardScaler()], ['clf', clf1]])
            pipe3 = Pipeline([['ms', MinMaxScaler()], ['clf', clf3]])
            all_clf = [pipe1, pipe3]
            clf_labels = ['Logistic Regression', 'SVM linear']
            cv = x_train.shape[0] if df_cv > x_train.shape[0] else df_cv
            for clf, label in zip(all_clf, clf_labels):
                cv_roc = cross_val_score(estimator=clf, X=x, y=ss1, cv=cv, scoring='roc_auc')
                clf.fit(x_train, y_train)
                y_pre_tr = clf.predict_proba(x_train)[:, 1]
                y_pre_ts = clf.predict_proba(x_test)[:, 1]
                fpr_tr, tpr_tr, thresholds_tr = roc_curve(y_true=y_train, y_score=y_pre_tr)
                fpr, tpr, thresholds = roc_curve(y_true=y_test, y_score=y_pre_ts)
                roc_auc_tr = auc(x=fpr_tr, y=tpr_tr)
                roc_auc = auc(x=fpr, y=tpr)
                if label == 'Logistic Regression':
                    l_cv_list.append(cv_roc.mean())
                    l_tr_list.append(roc_auc_tr)
                    l_ts_list.append(roc_auc)
                else:
                    ml_cv_list.append(cv_roc.mean())
                    ml_tr_list.append(roc_auc_tr)
                    ml_ts_list.append(roc_auc)
        if tn == 45:
            list_45l += l_ts_list
        if tn == 50:
            list_50l += l_ts_list
        print()
        df_s = pd.DataFrame()
        df_s['mir_number'] = pd.Series(fn_list)
        df_s['lr_cv10_roc'] = pd.Series(l_cv_list)
        df_s['lr_train_roc'] = pd.Series(l_tr_list)
        df_s['lr_test_roc'] = pd.Series(l_ts_list)
        df_s['svl_cv10_roc'] = pd.Series(ml_cv_list)
        df_s['svl_train_roc'] = pd.Series(ml_tr_list)
        df_s['svl_test_roc'] = pd.Series(ml_ts_list)
        lines = [':', '-.', '--']
        colors = ['red', 'blue', 'green']
        col_l = ['lr_cv10_roc', 'lr_train_roc', 'lr_test_roc']

        for ct, ls, clr in zip(col_l, lines, colors):
            ax[j, i].plot(df_s['mir_number'], df_s[ct], color=clr, linestyle=ls, label='Logistic Regression {}'.
                          format(ct.split('_', 1)[1]))
        ax[j, i].legend(loc='lower right')
        ax[j, i].set_title('Logistic Regression: Top[{}]'.format(tn))
        ax[j, i].grid()
        ax[j, i].set_xlabel('Feature numbers')
        ax[j, i].set_ylabel('ROC AUC')
        ax[j, i].set_yticks(np.arange(0.0, 1.05, 0.1))

        if tn in [45, 50, 55]:
            col_s = ['svl_cv10_roc', 'svl_train_roc', 'svl_test_roc']
            for ct, ls, clr in zip(col_s, lines, colors):
                ax[2, i].plot(df_s['mir_number'], df_s[ct], color=clr, linestyle=ls, label='Linear SVM {}'.
                              format(ct.split('_', 1)[1]))
                ax[2, i].legend(loc='lower right')
                ax[2, i].set_title('Linear SVM: Top[{}]'.format(tn))
                ax[2, i].grid()
                ax[2, i].set_xlabel('Feature numbers')
                ax[2, i].set_ylabel('ROC AUC')
                ax[2, i].set_yticks(np.arange(0.0, 1.05, 0.1))
        i += 1
        if i == 3:
            j += 1
            i = 0
    plt.tight_layout()
    plt.savefig(os.path.join(out_path, title_n))
    plt.close()
    return list_45l, list_50l


def mvc(df_xi, ss_y, title_n, out_path, ts, f_number, eps, df_cv):
    r_path = os.path.join(out_path, 'roc_curve')
    prepare_output_dir(r_path)
    pre_path = 'pos' if title_n.endswith('positive') else 'neg' if title_n.endswith('negative') else 'all'
    pre_path2 = os.path.join(r_path, pre_path)
    prepare_output_dir(pre_path2)
    df_x = df_xi.copy()
    fn_list, l_tr_list, l_cv_list, l_ts_list, m_tr_list, m_cv_list, m_ts_list = [], [], [], [], [], [], []
    ml_tr_list, ml_cv_list, ml_ts_list = [], [], []
    tmp1_list, tmp2_list = [], []
    print('Now running F_number:')
    for fn in range(1, f_number + 1):
        print(fn, end=',')
        fn_list.append(fn)
        x, feature_string, f_l = select_r(df_x, ss_y, fn, eps)
        tmp1_list.append('{},{},{}'.format(fn, f_l, feature_string))
        x_train, x_test, y_train, y_test = train_test_split(x, ss_y, test_size=ts, random_state=rs_stat)
        clf1 = LogisticRegression(penalty='l2', C=0.001, random_state=0)
        clf3 = SVC(kernel='linear', probability=True, random_state=0)
        pipe1 = Pipeline([['sc', StandardScaler()], ['clf', clf1]])
        pipe3 = Pipeline([['mc', MinMaxScaler()], ['clf', clf3]])
        all_clf = [pipe1, pipe3]
        clf_labels = ['Logistic Regression', 'SVM linear']
        cv = x_train.shape[0] if df_cv > x_train.shape[0] else df_cv
        for clf, label in zip(all_clf, clf_labels):
            cv_roc = cross_val_score(estimator=clf, X=x, y=ss_y, cv=cv, scoring='roc_auc')
            lf = title_n.split('_')
            clf.fit(x_train, y_train)
            y_pre_tr = clf.predict_proba(x_train)[:, 1]
            y_pre_ts = clf.predict_proba(x_test)[:, 1]
            fpr_tr, tpr_tr, thresholds_tr = roc_curve(y_true=y_train, y_score=y_pre_tr)
            fpr, tpr, thresholds = roc_curve(y_true=y_test, y_score=y_pre_ts)
            roc_auc_tr = auc(x=fpr_tr, y=tpr_tr)
            roc_auc = auc(x=fpr, y=tpr)
            n_tr, acc_tr, sen_tr, spe_tr, pre_tr = con_m(y_train, clf.predict(x_train))
            n_ts, acc_ts, sen_ts, spe_ts, pre_ts = con_m(y_test, clf.predict(x_test))
            tmp2_list.append('{},{},{},[{}],{:.2f} +/- {:.2f},{},{},{},{},{},{},{},{},{},{},{},{}'.
                             format(f_l, '_'.join(lf[0:3]), '_'.join(lf[3:]), label, cv_roc.mean(), cv_roc.std(), n_tr,
                                    roc_auc_tr, acc_tr, pre_tr, sen_tr, spe_tr, n_ts, roc_auc, acc_ts, pre_ts, sen_ts,
                                    spe_ts))
            if label == 'Logistic Regression':
                l_cv_list.append(cv_roc.mean())
                l_tr_list.append(roc_auc_tr)
                l_ts_list.append(roc_auc)
            else:
                ml_cv_list.append(cv_roc.mean())
                ml_tr_list.append(roc_auc_tr)
                ml_ts_list.append(roc_auc)
            if fn >= 10:
                nx = x.copy()
                ny = ss_y.values
                cvf = StratifiedKFold(n_splits=cv)
                tprs, aucs = [], []
                mean_fpr = np.linspace(0, 1, 100)
                plt.rcParams["figure.figsize"] = [12, 9]
                for train, test in cvf.split(nx, ny):
                    md = clf.fit(nx[train], ny[train])
                    pb = md.predict_proba(nx[test])
                    fpr_cv, tpr_cv, thresholds_cv = metrics.roc_curve(ny[test], pb[:, 1])
                    tprs.append(interp(mean_fpr, fpr_cv, tpr_cv))
                    tprs[-1][0] = 0.0
                    roc_auc_cv = metrics.auc(fpr_cv, tpr_cv)
                    aucs.append(roc_auc_cv)
                mean_tpr = np.mean(tprs, axis=0)
                mean_tpr[-1] = 1.0
                mean_auc = metrics.auc(mean_fpr, mean_tpr)
                std_auc = np.std(aucs)
                plt.plot(mean_fpr, mean_tpr, color='b', label=r'CV Mean ROC (AUC = %0.2f $\pm$ %0.2f)' %
                                                              (mean_auc, std_auc), lw=2, alpha=.8)
                std_tpr = np.std(tprs, axis=0)
                # tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
                # tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
                # plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2, label=r'$\pm$ 1 std. dev.')
                plt.plot(fpr, tpr, color='r', linestyle=':', label='{} Test (AUC = {:.2f})'.format(label, roc_auc),
                         lw=2, alpha=.8)
                plt.plot(fpr_tr, tpr_tr, color='g', linestyle='-.', label='{} Train (auc = {:.2f})'
                         .format(label, roc_auc_tr), lw=2, alpha=.8)
                plt.legend(loc='lower right')
                plt.plot([0, 1], [0, 1], linestyle='--', color='gray', linewidth=2)
                plt.xlim([-0.1, 1.1])
                plt.ylim([-0.1, 1.1])
                plt.title('Receiver operating characteristic of [{}]'.format(label))
                plt.grid()
                plt.xlabel('False Positive Rate (1 - Specificity)')
                plt.ylabel('True Positive Rate (Sensitivity)')
                f_name = title_n.replace(' ', '_').replace('(', '').replace(')', '')
                plt.tight_layout()
                plt.savefig(os.path.join(pre_path2, '{}_{}_{}'.format(f_name, label.replace(' ', '_'), fn)))
                plt.close()
    print()
    df_s = pd.DataFrame()
    df_s['mir_number'] = pd.Series(fn_list)
    df_s['lr_cv10_roc'] = pd.Series(l_cv_list)
    df_s['lr_train_roc'] = pd.Series(l_tr_list)
    df_s['lr_test_roc'] = pd.Series(l_ts_list)
    df_s['svl_cv10_roc'] = pd.Series(ml_cv_list)
    df_s['svl_train_roc'] = pd.Series(ml_tr_list)
    df_s['svl_test_roc'] = pd.Series(ml_ts_list)
    fn_check_plot(df_s, out_path, title_n)
    return tmp1_list, '\n'.join(tmp2_list)


def plot_tt(dfm_in, ssp_in, ssp, path_o, phe, suf):
    o_dir = os.path.join(path_o, 'split_plot')
    prepare_output_dir(o_dir)
    dfm2 = dfm_in.copy()
    ssp2 = ssp_in.copy()
    x_train, x_test, y_train, y_test = train_test_split(dfm2, ssp2, test_size=0.4, random_state=rs_stat)
    ss_train = ssp[y_train.index]
    ss_test = ssp[y_test.index]
    ss_train.sort_values(inplace=True)
    ss_test.sort_values(inplace=True)
    ap = ss_train.plot(color='g', linestyle=':', label='train')
    ss_test.plot(color='b', linestyle='-.', label='test')
    ap.axhline(y=ssp.mean(), color='r', linestyle='--', label='mean')
    # ap.axhline(y=11, color='r', linestyle='--', label='split line')
    plt.legend(loc='lower right')
    plt.tight_layout()
    plt.savefig(os.path.join(o_dir, '{}_{}_rs{}_1'.format(phe, suf, rs_stat)))
    plt.close()

    s_raw = pd.concat([ss_train, ss_test])
    s_raw.sort_values(inplace=True)
    s_t = pd.Series(0, index=ss_train.index)
    s_s = pd.Series(1, index=ss_test.index)
    s_binary = pd.concat([s_t, s_s])
    df1 = pd.DataFrame()
    df1['raw'] = s_raw
    df1['bin'] = s_binary
    df1['index'] = range(1, 1 + len(s_raw))
    df1.reset_index(inplace=True)
    df_train = df1[df1['bin'] == 0]
    df_test = df1[df1['bin'] == 1]
    ax = df_train.plot(kind='scatter', x='index', y='raw', color='b', label='train', title='Training/test distribution')
    df_test.plot(kind='scatter', x='index', y='raw', color='r', label='test', ax=ax)
    ax.axhline(y=ssp.mean(), color='r', linestyle='--', label='mean')
    # ax.axhline(y=11, color='r', linestyle='--', label='split line')
    plt.grid()
    plt.xlabel('rice sample index')
    plt.ylabel(phe)
    plt.tight_layout()
    plt.savefig(os.path.join(o_dir, '{}_{}_rs{}_2'.format(phe, suf, rs_stat)))
    plt.close()


def pn(df_mir, df_phe, r_p, phenotype, f_prefix, na_suffix, top_cn, top_n, c_th, fn):
    r_dir = os.path.join(r_p, '{}_{}_c{}t{}th{}'.format(f_prefix, na_suffix, top_cn, top_n, c_th))
    prepare_output_dir(r_dir)
    dfm = df_mir.copy()
    dfp = df_phe.copy()
    ssp = dfp[phenotype]
    print('na sum:{}_{}:'.format(f_prefix, na_suffix), sum(ssp.isnull()), dfm.shape, dfp.shape)
    ssp = ssp.dropna()
    phe_m = phenotype.replace('(', '').replace(')', '').replace(' ', '_')
    c_path = os.path.join(r_dir, '{}_corT{}.csv'.format(f_prefix, top_cn))
    if not os.path.exists(c_path):
        df_c, df_p = rice_corr(dfm, dfp, top_cn, 'pearson')
        extract_rc(c_path, df_c, locate(df_c, 'c', 0.01), 'correlation')
        extract_rc(c_path.replace('cor', 'pvl'), df_p, locate(df_c, 'p', 1), 'p_value')
    ss1 = s_top_gt(ssp, top_n, True)
    p2gp, p2gm, m2c = cor_dict_get(c_path, c_th)
    df_mp = dfm.loc[ss1.index, p2gp.get(phenotype)]
    df_mm = dfm.loc[ss1.index, p2gm.get(phenotype)]
    df_ma = pd.concat([df_mp, df_mm], 1)

    c_px = '{}_{}_c{}t{}'.format(f_prefix, na_suffix, top_cn, top_n)
    p_dir = os.path.join(r_dir, phe_m)
    prepare_output_dir(p_dir)
    plot_tt(df_mp, ss1, ssp, p_dir, phe_m, c_px)

    #mir_p, roc_p = mvc(df_mp, ss1, '{}_{}_positive'.format(c_px, phenotype), p_dir, 0.4, fn, 300, 10)
    mir_m, roc_m = mvc(df_mm, ss1, '{}_{}_negative'.format(c_px, phenotype), p_dir, 0.4, fn, 300, 10)
    mir_a, roc_a = mvc(df_ma, ss1, '{}_{}_all'.format(c_px, phenotype), p_dir, 0.4, fn, 300, 10)
    for m, r, p in zip([mir_m, mir_a], [roc_m, roc_a], ['neg', 'all']):
        with open(os.path.join(p_dir, '{}_{}_{}_mir.csv'.format(c_px, phe_m, p)), 'w') as o1, \
                open(os.path.join(p_dir, '{}_{}_{}_roc.csv'.format(c_px, phe_m, p)), 'w') as o2:
            o2.write('Feature_number,table_process,Phenotype_mir_type,Classifier,10-fold cv,Training samples,ROC_AUC,'
                     'Accuracy,Precision,Sensitivity,Specificity,Testing samples,ROC_AUC,Accuracy,Precision,'
                     'Sensitivity,Specificity\n')
            o1.write('\n'.join(m) + '\n')
            o2.write(r + '\n')


def feature_score4(dfx, ssy_n, f_string, m2c, path_o, title_n, f_number, phe, top_n):
    o_path = os.path.join(path_o, 'feature_score')
    prepare_output_dir(o_path)
    pre_path = 'pos' if title_n.endswith('pos') else 'neg' if title_n.endswith('neg') else 'all'
    pre_path2 = os.path.join(o_path, pre_path)
    prepare_output_dir(pre_path2)

    dfs = dfx[f_string.split(',')]
    mc_list = [float(m2c.get(x)) for x in dfs.columns]
    ssc = pd.Series(mc_list, index=dfs.columns)
    fn_list = ['def', 'rev', 'non']
    pv_list = []
    for fn in fn_list:
        ssc = ssc if fn == 'def' else ssc*(-1) if fn == 'rev' else 1

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
                      .format(title_n, f_number), ax=ax[0])
        ax[0].axvline(x=(top_n * 2) / 3 + 0.5, color='g', linestyle=':')
        ax[0].axvline(x=(top_n * 4) / 3 + 0.5, color='g', linestyle=':')
        ax[0].text((top_n * 2) / 3 + 0.5, 0.5 * (df1['scores'].min() + df1['scores'].max()), '1/3 score', rotation=90)
        ax[0].text((top_n * 4) / 3 + 0.5, 0.5 * (df1['scores'].min() + df1['scores'].max()), '2/3 score', rotation=90)
        ax[0].set_ylabel('feature scores', color='b')
        ax[0].tick_params('y')
        ax[0].grid()
        df1.plot(kind='scatter', x='index', y=phe, color='r', label='phenotype', ax=ax[1])
        df1.plot.line(x='index', y='trend', ax=ax[1])
        ax[1].axvline(x=(top_n * 2) / 3 + 0.5, color='g', linestyle=':')
        ax[1].axvline(x=(top_n * 4) / 3 + 0.5, color='g', linestyle=':')
        ax[1].text((top_n * 2) / 3 + 0.5, 0.5 * (df1[phe].min() + df1[phe].max()), 'demarcation 1', rotation=90)
        ax[1].text((top_n * 4) / 3 + 0.5, 0.5 * (df1[phe].min() + df1[phe].max()), 'demarcation 2', rotation=90)
        ax[1].set_ylabel(phe, color='r')
        ax[1].tick_params('y')
        ax[1].grid()
        plt.tight_layout()
        plt.savefig(os.path.join(pre_path2, '{}_{}_{}'.format(title_n, f_number, fn)))
        plt.close()
        c1 = (top_n * 2) // 3
        c2 = (top_n * 4) // 3
        ss1 = df1[phe][:c1]
        ss2 = df1[phe][c1:c2]
        ss3 = df1[phe][c2:]
        tt1, pp1 = scipy_ttest_ind(ss1, ss2, False)
        tt2, pp2 = scipy_ttest_ind(ss2, ss3, False)
        tt3, pp3 = scipy_ttest_ind(ss1, ss3, False)
        pv_list.append('{}p1p2,{},{},{}'.format(fn, ss1.mean(), ss2.mean(), pp1))
        pv_list.append('{}p2p3,{},{},{}'.format(fn, ss2.mean(), ss3.mean(), pp2))
        pv_list.append('{}p1p3,{},{},{}'.format(fn, ss1.mean(), ss3.mean(), pp3))
    with open(os.path.join(pre_path2, '{}_{}.csv'.format(title_n, f_number)), 'w') as out_f:
        out_f.write('type,subset1,subset2,p_value\n')
        out_f.write('\n'.join(pv_list) + '\n')


def feature_score3(dfx, ssy_n, f_string, m2c, path_o, title_n, f_number, phe, top_n):
    o_path = os.path.join(path_o, 'feature_score')
    prepare_output_dir(o_path)
    pre_path = 'pos' if title_n.endswith('pos') else 'neg' if title_n.endswith('neg') else 'all'
    pre_path2 = os.path.join(o_path, pre_path)
    prepare_output_dir(pre_path2)

    dfs = dfx[f_string.split(',')]
    mc_list = [float(m2c.get(x)) for x in dfs.columns]
    ssc = pd.Series(mc_list, index=dfs.columns)
    fn_list = ['def', 'rev', 'non']
    pv_list = []
    for fn in fn_list:
        ssc = ssc if fn == 'def' else ssc*(-1) if fn == 'rev' else 1

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
                      .format(title_n, f_number), ax=ax[0])
        ax[0].axvline(x=(top_n * 2 + 1)/2, color='g', linestyle=':')
        ax[0].text((top_n * 2 + 1)/2, 0.5 * (df1['scores'].min() + df1['scores'].max()), 'median score', rotation=90)
        ax[0].set_ylabel('feature scores', color='b')
        ax[0].tick_params('y')
        ax[0].grid()
        df1.plot(kind='scatter', x='index', y=phe, color='r', label='phenotype', ax=ax[1])
        df1.plot.line(x='index', y='trend', ax=ax[1])
        ax[1].axvline(x=(top_n * 2 + 1)/2, color='g', linestyle=':')
        ax[1].text((top_n * 2 + 1)/2, 0.5*(df1[phe].min() + df1[phe].max()), 'median score', rotation=90)
        ax[1].set_ylabel(phe, color='r')
        ax[1].tick_params('y')
        ax[1].grid()
        plt.tight_layout()
        plt.savefig(os.path.join(pre_path2, '{}_{}_{}'.format(title_n, f_number, fn)))
        plt.close()
        # if fn == 'non':
        #     df1.to_csv(os.path.join(pre_path2, '{}_{}_{}z.csv'.format(title_n, f_number, fn)))
        ss1 = df1[phe][:top_n]
        ss2 = df1[phe][-top_n:]
        tt, pp = scipy_ttest_ind(ss1, ss2, False)
        pv_list.append('{},{},{},{},{}'.format(fn, ss1.mean(), ss2.mean(), pp, tt))
    with open(os.path.join(pre_path2, '{}_{}.csv'.format(title_n, f_number)), 'w') as out_f:
        out_f.write('type,subset1,subset2,p_value,t_statistic\n')
        out_f.write('\n'.join(pv_list) + '\n')


def out_feature_plot(dfx, ssy, ssy_n, f_string, m2c, m2p, path_o, title_n, f_number):
    o_path = os.path.join(path_o, 'Features')
    prepare_output_dir(o_path)
    pre_path = 'pos' if title_n.endswith('pos') else 'neg' if title_n.endswith('neg') else 'all'
    pre_path2 = os.path.join(o_path, pre_path)
    prepare_output_dir(pre_path2)
    lr_coe, lr_int, lr_score = linear_rg(dfx, ssy_n, f_string)
    df_p = t_test(dfx, ssy, f_string)
    df_p.to_csv(os.path.join(pre_path2, '{}_{}_pvl.csv'.format(title_n, f_number)))
    high_df = dfx[ssy == 1]
    low_df = dfx[ssy == 0]
    f_list = f_string.split(',')
    fig, axes = plt.subplots(math.ceil(f_number / 2), 2, figsize=(10, 2 * math.ceil(f_number / 2)))
    ax = axes.ravel()
    with open(os.path.join(pre_path2, '{}_{}_cor.csv'.format(title_n, f_number)), 'w') as out_c:
        i = 0
        out_c.write('miRNA,correlation,p_value,lr_coe(intercept {:.3f})\n'.format(lr_int))
        for m, c in zip(f_list, lr_coe):
            out_c.write('{},{},{},{}\n'.format(m, m2c.get(m, 'No_result'), m2p.get(m, 'No_result'), c))
            _, bins = np.histogram(dfx[m], bins=50)
            ax[i].hist(high_df[m], bins=bins, color='r', alpha=.5)
            ax[i].hist(low_df[m], bins=bins, color='b', alpha=.5)
            ax[i].set_title(m)
            ax[i].set_yticks(())
            i += 1
        ax[1].set_xlabel('Feature magnitude')
        ax[1].set_ylabel('Frequency')
        ax[1].legend(['High', 'Low'], loc='best')
        fig.tight_layout()
        plt.suptitle(title_n)
        plt.savefig(os.path.join(pre_path2, '{}_{}_mir'.format(title_n, f_number)))
        plt.close()
        out_c.write(',,lr_score,{}\n'.format(lr_score))


def linear_rg(dfx, ssy, f_string):
    df = dfx.copy()
    lr = LinearRegression()
    x = df[f_string.split(',')]
    lr.fit(x, ssy)
    return lr.coef_, lr.intercept_, lr.score(x, ssy)


def scipy_ttest_ind(s1, s2, var=False):
    return st.ttest_ind(s1, s2, equal_var=var)


def t_test(dfx, ssy, f_string):
    x = dfx[f_string.split(',')]
    y = ssy.copy()
    x0 = x.loc[y[y == 0].index, :]
    x1 = x.loc[y[y == 1].index, :]
    df1 = pd.DataFrame()
    df1['Low_mean'] = x0.mean()
    df1['High_mean'] = x1.mean()
    df1['FC(High/Low)'] = x1.mean() / x0.mean()
    p_list = []
    for m in df1.index:
        pv = scipy_ttest_ind(x0[m], x1[m])[1]
        p_list.append('{:.3g}'.format(pv))
    df1['ind_t_test'] = np.array(p_list)
    df1['Low_std'] = x0.std()
    df1['High_std'] = x1.std()
    return df1


def plot_f50_nap(n_l, a_l, p_l, path_o, tp):
    df1 = pd.DataFrame()
    df1['mir_number'] = pd.Series(list(range(1, 50)))
    df1['negative'] = pd.Series(n_l)
    df1['all'] = pd.Series(a_l)
    df1['positive'] = pd.Series(p_l)
    lines = ['-.', ':', '--']
    colors = ['blue', 'red', 'green']
    col_l = ['negative', 'all', 'positive']
    plt.rcParams["figure.figsize"] = [8, 6]
    for ct, ls, clr in zip(col_l, lines, colors):
        plt.plot(df1['mir_number'], df1[ct], color=clr, linestyle=ls, label='{} ROC AUC'.format(ct))
    plt.legend(loc='lower right')
    plt.title('Logistic regression of all, positive, and negative correlation miRNAs')
    plt.grid()
    plt.xticks(np.arange(0, max(df1['mir_number']) + 1, 5))
    plt.yticks(np.arange(0.0, 1.05, 0.1))
    plt.xlabel('Feature numbers')
    plt.tight_layout()
    plt.savefig(os.path.join(path_o, 'miRNA_comparison_{}'.format(tp)))
    plt.close()


def main():
    check1 = True
    test1 = True
    check2 = True
    check3 = True
    main_dir = './'
    dfr, n_p = df_open(os.path.join(main_dir, 'wm_all_q.csv'))
    dfr = dfr.drop(['type (H)', 'waxy (H)'], axis=1)
    dfr2 = dfr.dropna()
    dfm = dfr.iloc[:, :924]
    dfp = dfr.iloc[:, 924:]
    dfm2 = dfr2.iloc[:, :924]
    dfp2 = dfr2.iloc[:, 924:]
    ct = 50
    tn = 45
    th = 0.15
    fn = 50
    if check1:
        print('Step 1')
        c_list = ['Panicle Number (I)'] if test1 else dfp2.columns
        for phenotype_r in c_list:
            print('Phenotype: {}'.format(phenotype_r))
            #  pn(dfm, dfp, main_dir, phenotype_r, 'wmq', 'no_drop', ct, tn, th, fn)
            pn(dfm2, dfp2, main_dir, phenotype_r, 'wmq', 'drop', ct, tn, th, fn)
    if check2:
        print('Step 2')
        # d_list = ['drop', 'no_drop']
        d_list = ['drop']
        c_list = ['Panicle Number (I)'] if test1 else dfp2.columns
        for n in d_list:
            r_dir = os.path.join(main_dir, '{}_{}_c{}t{}th{}'.format('wmq', n, ct, tn, th))
            c_path = os.path.join(r_dir, 'wmq_corT{}.csv'.format(ct))
            p2gp, p2gm, m2c = cor_dict_get(c_path, th)
            m2p = p_dict_get(os.path.join(r_dir, 'wmq_pvlT{}.csv'.format(ct)), 0.8)
            for phenotype in c_list:
                print('Phenotype: {}'.format(phenotype))
                phe_m = phenotype.replace('(', '').replace(')', '').replace(' ', '_')
                p_dir = os.path.join(r_dir, phe_m)
                prepare_output_dir(p_dir)
                c_px = '{}_{}_c{}t{}'.format('wmq', n, ct, tn)
                ssp = dfp2[phenotype] if n == 'drop' else dfp[phenotype]
                ssp = ssp.dropna()
                ss1 = s_top_gt(ssp, tn, True)
                ssp = ssp[ss1.index]
                dfx = dfm2.loc[ss1.index, :] if n == 'drop' else dfm.loc[ss1.index, :]
                for c in ['neg', 'all']:
                    with open(os.path.join(p_dir, '{}_{}_{}_mir.csv'.format(c_px, phe_m, c))) as in_f:
                        for line in in_f:
                            lf = line.rstrip().split(',', maxsplit=2)
                            if 10 <= int(lf[1]) <= 30:
                                out_feature_plot(dfx, ss1, ssp, lf[2], m2c.get(phenotype), m2p.get(phenotype), p_dir,
                                                 '{}_{}_{}'.format(c_px, phe_m, c), int(lf[1]))
                                feature_score3(dfx, ssp, lf[2], m2c.get(phenotype), p_dir, '{}_{}_{}'.
                                               format(c_px, phe_m, c), int(lf[1]), phe_m, tn)
    if check3:
        print('Step 3')
        o_dir = './n_test'
        prepare_output_dir(o_dir)
        c_path = os.path.join(o_dir, 'wmq_corT50.csv')
        if not os.path.exists(c_path):
            df_c, df_p = rice_corr(dfm2, dfp2, 50, 'pearson')
            extract_rc(c_path, df_c, locate(df_c, 'c', 0.01), 'correlation')
        p2gp, p2gm, m2c = cor_dict_get(c_path, 0.15)
        ssp = dfp2['Panicle Number (I)']
        df_mp = dfm2.loc[:, p2gp.get('Panicle Number (I)')]
        df_mm = dfm2.loc[:, p2gm.get('Panicle Number (I)')]
        df_ma = pd.concat([df_mp, df_mm], 1)
        neg_45l, neg_50l = top_n_test(df_mm, ssp, 'Panicle_Number_I_negative', o_dir, 0.4, fn, 300, 10)
        all_45l, all_50l = top_n_test(df_ma, ssp, 'Panicle_Number_I_all', o_dir, 0.4, fn, 300, 10)
        pos_45l, pos_50l = top_n_test(df_mp, ssp, 'Panicle_Number_I_positive', o_dir, 0.4, fn, 300, 10)
        plot_f50_nap(neg_45l, all_45l, pos_45l, o_dir, 45)
        plot_f50_nap(neg_50l, all_50l, pos_50l, o_dir, 50)


if __name__ == '__main__':
    main()
