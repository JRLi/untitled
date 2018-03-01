#!/usr/bin/env python3.6
import pandas as pd
import numpy as np
import scipy.stats as st
import os
import math
import matplotlib.pyplot as plt
from collections import defaultdict
from sklearn.linear_model import LinearRegression
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
    return ss


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
    l_path = os.path.join(path_o, 'roc_lg')
    m_path = os.path.join(path_o, 'roc_mv')
    prepare_output_dir(l_path)
    prepare_output_dir(m_path)
    df1 = df_in.copy()
    for t, o_dir in zip(['Logistic Regression', 'Majority Voting'], [l_path, m_path]):
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
        plt.savefig(os.path.join(o_dir, '{}_{}'.format(f_name, t.replace(' ', '_'))))
        plt.gcf().clear()


def mvc(df_xi, ss_y, title_n, out_path, c_dict, ts, f_number, eps, df_cv):
    df_x = df_xi.copy()
    fn_list, l_cvc_list, l_cvu_list, l_roc_list, m_cvc_list, m_cvu_list, m_roc_list = [], [], [], [], [], [], []
    tmp1_list, tmp2_list = [], []
    print('Now running F_number:')
    for fn in range(1, f_number + 1):
        print(fn, end=',')
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
            lf = title_n.split('_')
            y_pre = clf.fit(x_train, y_train).predict_proba(x_test)[:, 1]
            fpr, tpr, thresholds = roc_curve(y_true=y_test, y_score=y_pre)
            roc_auc = auc(x=fpr, y=tpr)
            tmp2_list.append('{}\t{}\t{}\t[{}]\t{:.2f}\t+/- {:.2f}\t{:.2f}\t+/- {:.2f}\t{:.2f}'.
                             format(f_l, '_'.join(lf[0:3]), '_'.join(lf[3:]), label, cv_roc.mean(), cv_roc.std(),
                                    cv_acc.mean(), cv_acc.std(), roc_auc))
            if label == 'Logistic Regression':
                l_cvc_list.append(cv_acc.mean())
                l_cvu_list.append(cv_roc.mean())
                l_roc_list.append(roc_auc)
            else:
                m_cvc_list.append(cv_acc.mean())
                m_cvu_list.append(cv_roc.mean())
                m_roc_list.append(roc_auc)
    print()
    df_s = pd.DataFrame()
    df_s['mir_number'] = pd.Series(fn_list)
    df_s['lr_cv10_acc'] = pd.Series(l_cvc_list)
    df_s['lr_cv10_roc'] = pd.Series(l_cvu_list)
    df_s['lr_test_roc'] = pd.Series(l_roc_list)
    df_s['mv_cv10_acc'] = pd.Series(m_cvc_list)
    df_s['mv_cv10_roc'] = pd.Series(m_cvu_list)
    df_s['mv_test_roc'] = pd.Series(m_roc_list)
    fn_check_plot(df_s, out_path, title_n)
    return tmp1_list, '\n'.join(tmp2_list)


def pn_summary(ssp2, path_o, suf):
    ssp = ssp2.copy()
    with open(os.path.join(path_o, 'top_summary_{}'.format(suf)), 'w') as of1:
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


def plot_tt(dfm_in, ssp_in, ssp, path_o, phe, suf, rs=0):
    o_dir = os.path.join(path_o, 'split_plot')
    prepare_output_dir(o_dir)
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
    plt.savefig(os.path.join(o_dir, '{}_{}_rs{}_1'.format(phe, suf, rs)))
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
    plt.ylabel(phe)
    plt.savefig(os.path.join(o_dir, '{}_{}_rs{}_2'.format(phe, suf, rs)))
    plt.gcf().clear()


def pn(df_mir, df_phe, r_p, phenotype, f_prefix, na_suffix, top_cn, top_n, c_th):
    r_dir = os.path.join(r_p, '{}_{}_c{}t{}th{}'.format(f_prefix, na_suffix, top_cn, top_n, c_th))
    prepare_output_dir(r_dir)
    dfm = df_mir.copy()
    dfp = df_phe.copy()
    ssp = dfp[phenotype]
    print('na sum:{}_{}:'.format(f_prefix, na_suffix), sum(ssp.isnull()), dfm.shape, dfp.shape)
    ssp = ssp.dropna()
    phe_m = phenotype.replace('(', '').replace(')', '').replace(' ', '_')

    pn_summary(ssp, r_dir, '{}_{}_{}'.format(phe_m, f_prefix, na_suffix))

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
    print(df_mp.shape, df_mm.shape, df_ma.shape)
    # check top 50 sample split into training and test
    c_px = '{}_{}_c{}t{}'.format(f_prefix, na_suffix, top_cn, top_n)

    plot_tt(df_mp, ss1, ssp, r_dir, phe_m, c_px, 1)

    mir_p, roc_p = mvc(df_mp, ss1, '{}_{}_positive'.format(c_px, phenotype), r_dir, m2c.get(phenotype), 0.4, 30, 300, 10)
    mir_m, roc_m = mvc(df_mm, ss1, '{}_{}_negative'.format(c_px, phenotype), r_dir, m2c.get(phenotype), 0.4, 30, 300, 10)
    mir_a, roc_a = mvc(df_ma, ss1, '{}_{}_p_and_n'.format(c_px, phenotype), r_dir, m2c.get(phenotype), 0.4, 30, 300, 10)
    for m, r, p in zip([mir_p, mir_m, mir_a], [roc_p, roc_m, roc_a], ['pos', 'neg', 'all']):
        with open(os.path.join(r_dir, '{}_{}_{}_mir'.format(c_px, phe_m, p)), 'w') as o1, \
                open(os.path.join(r_dir, '{}_{}_{}_roc'.format(c_px, phe_m, p)), 'w') as o2:
            o1.write('\n'.join(m) + '\n')
            o2.write(r + '\n')


def out_feature_plot(dfx, ssy, f_string, m2c, path_o, title_n, f_number):
    f_name = title_n.replace(' ', '_').replace('(', '').replace(')', '')
    o_path = os.path.join(path_o, 'Features')
    prepare_output_dir(o_path)
    high_df = dfx[ssy == 1]
    low_df = dfx[ssy == 0]
    f_list = f_string.split(',')
    fig, axes = plt.subplots(math.ceil(f_number / 2), 2, figsize=(10, 2 * math.ceil(f_number / 2)))
    ax = axes.ravel()
    with open(os.path.join(o_path, '{}_{}_cor.csv'.format(f_name, f_number)), 'w') as out_c:
        i = 0
        for m in f_list:
            out_c.write('{},{}\n'.format(m, m2c.get(m, 'No_result')))
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
        plt.savefig(os.path.join(o_path, '{}_{}_mir'.format(f_name, f_number)))
        #plt.gcf().clear()
        plt.close()


def linear_rg(dfx, ssy, f_string, path_o, title_n, f_number):
    pass


def main():
    x = False
    main_dir = 'E:/StringTemp/Project_Rice'
    phenotype_list = ['Panicle Number (I)']
    dfr, n_p = df_open(os.path.join(main_dir, 'wm_df129_q.csv'))
    dfr = dfr.drop(['type (H)', 'waxy (H)'], axis=1)  # df2, drop all na first
    dfr2 = dfr.dropna()
    dfm = dfr.iloc[:, :924]
    dfp = dfr.iloc[:, 924:]
    dfm2 = dfr2.iloc[:, :924]
    dfp2 = dfr2.iloc[:, 924:]
    ct = 50
    tn = 45
    th = 0.15
    if x:
        for phenotype_r in phenotype_list:
            pn(dfm, dfp, main_dir, phenotype_r, 'wmq', 'no_drop', ct, tn, th)
            pn(dfm2, dfp2, main_dir, phenotype_r, 'wmq', 'drop', ct, tn, th)
    else:
        for n in ['drop', 'no_drop']:
            r_dir = os.path.join(main_dir, '{}_{}_c{}t{}th{}'.format('wmq', n, ct, tn, th))
            c_path = os.path.join(r_dir, 'wmq_corT{}.csv'.format(ct))
            p2gp, p2gm, m2c = cor_dict_get(c_path, th)
            for phenotype in phenotype_list:
                phe_m = phenotype.replace('(', '').replace(')', '').replace(' ', '_')
                c_px = '{}_{}_c{}t{}'.format('wmq', n, ct, tn)
                ssp = dfp2[phenotype] if n == 'drop' else dfp[phenotype]
                ssp = ssp.dropna()
                ss1 = s_top_gt(ssp, tn, True)
                dfx = dfm2.loc[ss1.index, :] if n == 'drop' else dfm.loc[ss1.index, :]
                for c in ['pos', 'neg', 'all']:
                    with open(os.path.join(r_dir, '{}_{}_{}_mir'.format(c_px, phe_m, c))) as in_f:
                        for line in in_f:
                            lf = line.rstrip().split('\t')
                            if int(lf[1]) <= 20:
                                out_feature_plot(dfx, ss1, lf[2], m2c.get(phenotype), r_dir,
                                                 '{}_{}_{}'.format(c_px, phe_m, c), int(lf[1]))


if __name__ == '__main__':
    main()
