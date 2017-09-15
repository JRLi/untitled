#!/usr/bin/env python3
import os
import pandas as pd

dir_loc = 'D:/Project/circ_bicluster/InCoB_2017_attachment_71/InCoB article 71/supplemental_and_cover_letter/Additional_File_2'


def open_df(in_path, direct='n'):
    fpath, fname = os.path.split(in_path)
    fbase, fext = os.path.splitext(fname)
    df = pd.read_csv(in_path, index_col=0) if fext == '.csv' else pd.read_table(in_path, index_col=0)
    #print('{}: Transpose: {}'.format(fbase, direct))
    if direct == 't':
        df = df.transpose()
    return df, fbase


def main():
    file_list = os.listdir(dir_loc)
    print(len(file_list))
    m_list, g_list, s_list= [], [], []
    m2s_dict, g2s_dict, s2s_dict = {}, {}, {}
    for file in file_list:
        df1, df1_b = open_df(os.path.join(dir_loc, file))
        print(df1.shape[0], df1.shape[1], df1.shape[0] * df1.shape[1])
        m_list.append(df1.shape[0] * df1.shape[1])
        g_list.append(df1.shape[0])
        s_list.append(df1.shape[1])
        m2s_dict[df1.shape[0] * df1.shape[1]] = df1.shape
        g2s_dict[df1.shape[0]] = df1.shape
        s2s_dict[df1.shape[1]] = df1.shape
    m_list.sort()
    g_list.sort()
    s_list.sort()
    print(g2s_dict.get(g_list[0]), g2s_dict.get(g_list[-1]))
    print(s2s_dict.get(s_list[0]), s2s_dict.get(s_list[-1]))
    print(m2s_dict.get(m_list[0]), m2s_dict.get(m_list[-1]))

if __name__ == '__main__':
    main()