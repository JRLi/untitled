#!/usr/bin/env python3.6
import pandas as pd
import numpy as np
import scipy.stats as st
import os
import matplotlib.pyplot as plt
from collections import defaultdict
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFE
from sklearn.neighbors import KNeighborsClassifier
from sklearn.pipeline import Pipeline
from sklearn.tree import DecisionTreeClassifier
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_curve
from sklearn.metrics import auc
from sklearn.ensemble import VotingClassifier


def df_open(path_in, direct='n'):
    f_p, f_n = os.path.split(path_in)
    n_p, n_s = os.path.splitext(f_n)
    df1 = pd.read_csv(path_in, index_col=0) if n_s == '.csv' else pd.read_table(path_in, index_col=0)
    if direct != 'n':
        df1 = df1.transpose()
    return df1, f_p, n_p


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
    return ss


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


def fn_check_plot(df_in, path_o, title_n):
    df1 = df_in.copy()
    print(df1)
    for t in ['Logistic Regression', 'Majority Voting']:
        lines = [':', '-.', '--']
        colors = ['green', 'blue', 'orange']
        col_l = ['lr_cv10_acc', 'lr_cv10_roc', 'lr_test_roc'] if t == 'Logistic Regression' \
            else ['mv_cv10_acc', 'mv_cv10_roc', 'mv_test_roc']
        for ct, ls, clr in zip(col_l, lines, colors):
            plt.plot(df1['mir_number'], df1[ct], color=clr, linestyle=ls, label=ct.replace(ct[:2], t))
        plt.legend(loc='lower right')
        plt.title(title_n)
        plt.grid()
        plt.xlabel('Feature numbers')
        f_name = title_n.replace(' ', '_').replace('(', '').replace(')', '')
        plt.savefig(os.path.join(path_o, '{}_{}'.format(f_name, t.replace(' ', '_'))))
        plt.gcf().clear()


def mvc(df_xi, ss_y, title_n, out_path, ts, f_number, eps, df_cv):
    df_x = df_xi.copy()
    fn_list, l_cvc_list, l_cvu_list, l_roc_list, m_cvc_list, m_cvu_list, m_roc_list = [], [], [], [], [], [], []
    tmp1_list, tmp2_list = [], []
    for fn in range(1, f_number + 1):
        print('Now running F_number: {}'.format(fn))
        fn_list.append(fn)
        x, feature_string, f_l = select_r(df_x, ss_y, fn, eps)
        tmp1_list.append('{}\t{}\t{}'.format(fn, f_l, feature_string))
        x_train, x_test, y_train, y_test = train_test_split(x, ss_y, test_size=ts, random_state=1)
        clf1 = LogisticRegression(penalty='l2', C=0.001, random_state=0)
        clf2 = DecisionTreeClassifier(max_depth=3, criterion='entropy', random_state=0)
        clf3 = KNeighborsClassifier(n_neighbors=5, p=2, metric='minkowski')
        pipe1 = Pipeline([['sc', StandardScaler()], ['clf', clf1]])
        pipe3 = Pipeline([['sc', StandardScaler()], ['clf', clf3]])
        mv_clf = VotingClassifier([('lr', pipe1), ('dt', clf2), ('knn', pipe3)], voting='soft')
        all_clf = [pipe1, mv_clf]
        clf_labels = ['Logistic Regression', 'Majority Voting']
        cv = x_train.shape[0] if df_cv > x_train.shape[0] else df_cv
        for clf, label in zip(all_clf, clf_labels):
            cv_roc = cross_val_score(estimator=clf, X=x, y=ss_y, cv=cv, scoring='roc_auc')
            cv_acc = cross_val_score(estimator=clf, X=x, y=ss_y, cv=cv)
            lf = title_n.split('_', 1)
            y_pre = clf.fit(x_train, y_train).predict_proba(x_test)[:, 1]
            fpr, tpr, thresholds = roc_curve(y_true=y_test, y_score=y_pre)
            roc_auc = auc(x=fpr, y=tpr)
            tmp2_list.append('{}\t{}\t{}\t[{}]\t{:.2f}\t+/- {:.2f}\t{:.2f}\t+/- {:.2f}\t{:.2f}'.
                             format(f_l, lf[0], lf[1], label, cv_roc.mean(), cv_roc.std(), cv_acc.mean(), cv_acc.std(),
                                    roc_auc))
            if label == 'Logistic Regression':
                l_cvc_list.append(cv_acc.mean())
                l_cvu_list.append(cv_roc.mean())
                l_roc_list.append(roc_auc)
            else:
                m_cvc_list.append(cv_acc.mean())
                m_cvu_list.append(cv_roc.mean())
                m_roc_list.append(roc_auc)
    df_s = pd.DataFrame()
    df_s['mir_number'] = pd.Series(fn_list)
    df_s['lr_cv10_acc'] = pd.Series(l_cvc_list)
    df_s['lr_cv10_roc'] = pd.Series(l_cvu_list)
    df_s['lr_test_roc'] = pd.Series(l_roc_list)
    df_s['mv_cv10_acc'] = pd.Series(m_cvc_list)
    df_s['mv_cv10_roc'] = pd.Series(m_cvu_list)
    df_s['mv_test_roc'] = pd.Series(m_roc_list)
    fn_check_plot(df_s, out_path, title_n)
    return '\n'.join(tmp1_list), '\n'.join(tmp2_list)


def pn_summary(ssp2, path_o):
    ssp = ssp2.copy()
    with open(os.path.join(path_o, 'top_summary_for_PNI'), 'w') as of1:
        of1.write('TopN\tLowN\tHighN\tLowRange\tHighRange\n')
        for tn in [0, 60, 55, 50, 45, 40]:
            ssp_t = s_top_gt(ssp, tn, True)
            of1.write('{}\t{}\t{}\t'.format(tn, len(ssp_t[ssp_t == 0]), len(ssp_t[ssp_t == 1])))
            ssp3 = ssp[ssp_t.index]
            ss_l = ssp3[ssp_t == 0]
            ss_h = ssp3[ssp_t == 1]
            ss_l.sort_values(inplace=True)
            ss_h.sort_values(inplace=True)
            of1.write('{} ~ {}\t{} ~ {}\n'.format(ss_l[0], ss_l[-1], ss_h[0], ss_h[-1]))


def plot_tt(dfm_in, ssp_in, ssp, path_o, rs):
    dfm2 = dfm_in.copy()
    ssp2 = ssp_in.copy()
    x_train, x_test, y_train, y_test = train_test_split(dfm2, ssp2, test_size=0.4, random_state=rs)
    ss_train = ssp[y_train.index]
    ss_test = ssp[y_test.index]
    ss_train.sort_values(inplace=True)
    ss_test.sort_values(inplace=True)
    ap = ss_train.plot(color='g', linestyle=':', label='train')
    ss_test.plot(color='b', linestyle='-.', label='test')
    ap.axhline(y=ssp.mean(), color='r', linestyle='--', label='mean')
    plt.legend(loc='lower right')
    plt.savefig(os.path.join(path_o, 'Panicle_number_rs{}_1'.format(rs)))
    plt.gcf().clear()
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
    plt.grid()
    plt.xlabel('rice sample index')
    plt.ylabel('Panicle number')
    plt.savefig(os.path.join(path_o, 'Panicle_number_rs{}_2'.format(rs)))
    plt.gcf().clear()


def pn(df_mir, df_phe, r_p):
    dfm = df_mir.copy()
    dfp = df_phe.copy()
    ssp = dfp['Panicle Number (I)']
    print(sum(ssp.isnull()))
    ssp = ssp.dropna()
    dfm2 = dfm.loc[ssp.index, :]
    pn_summary(ssp, r_p)
    c_path = 'E:/StringTemp/Project_Rice/wm_df129_q_corT{}.csv'.format(50)
    if not os.path.exists(c_path):
        df_c, df_p = rice_corr(dfm, dfp, 50, 'pearson')
        extract_rc(c_path, df_c, locate(df_c, 'c', 0.01), 'correlation')
        extract_rc(c_path.replace('cor', 'pvl'), df_p, locate(df_c, 'p', 1), 'p_value')
    p2gp, p2gm = cor_dict_get(c_path, 0.1)
    ss1 = s_top_gt(ssp, 45, True)
    df_mp = dfm2.loc[ss1.index, p2gp.get('Panicle Number (I)')]
    df_mm = dfm2.loc[ss1.index, p2gm.get('Panicle Number (I)')]
    df_ma = pd.concat([df_mp, df_mm], 1)
    print(df_mp.shape, df_mm.shape, df_ma.shape)
    plot_tt(df_mp, ss1, ssp, r_p, 1)   # check top 50 sample split into training and test
    #mir_p, roc_p = mvc(df_mp, ss1, '{}_positive'.format('Panicle Number (I)'), r_p, 0.4, 30, 300, 10)
    mir_m, roc_m = mvc(df_mm, ss1, '{}_negative'.format('Panicle Number (I)'), r_p, 0.4, 30, 300, 10)
    with open(os.path.join(r_p, 'pn_mir'), 'w') as o1, open(os.path.join(r_p, 'pn_roc'), 'w') as o2:
        o1.write(mir_m)
        o2.write(roc_m)


def main():
    dfr, r_p, n_p = df_open('E:/StringTemp/Project_Rice/wm_df129_q.csv')
    dfm = dfr.iloc[:, :924]
    dfp = dfr.iloc[:, 924:]
    dfp = dfp.drop(['type (H)', 'waxy (H)'], axis=1)
    pn(dfm, dfp, r_p)


if __name__ == '__main__':
    main()
