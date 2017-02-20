#!/usr/bin/env python3.5
import numpy as np
import pandas as pd
import sys
filePath = 'E:/StringTemp/GDS/'
fileName = 'GDS3876.matrix'


def quantileNormalize(df_input):
    df = df_input.copy()
    #compute rank
    dic = {}
    for col in df:
        dic.update({col : sorted(df[col])})
    sorted_df = pd.DataFrame(dic)
    rank = sorted_df.mean(axis = 1).tolist()
    #sort
    for col in df:
        t = np.searchsorted(np.sort(df[col]), df[col])
        df[col] = [rank[i] for i in t]
    return df


def percentile(df_input, per_th = 10):
    df = df_input.copy()
    for i in range(0, len(df.columns)):
        qv = np.percentile(df.iloc[:, i], per_th)
        print(qv)
        df.iloc[:, i][df.iloc[:, i] < qv] = qv
        # print(df.iloc[:, i])
    return df


def df_mean_index(df_input):
    # gp1 = df1.groupby(df1.index)
    # df1 = gp1.mean()
    gp = df_input.groupby(df_input.index)
    return gp.mean()


def df_mean_columns(df_input):
    gp = df_input.groupby(df_input.columns, axis=1)
    return gp.mean()


def hgnc_dict(hgnc_path):
    h_dict = {}
    with open(hgnc_path, 'r', encoding='utf-8') as input:
        for line in input:
            stringField = line.split('\t', maxsplit=2)
            h_dict[stringField[1]] = stringField[0]
    return  h_dict


def main(argv=None):
    if argv is None:
        df1 = pd.read_table(filePath + 'GDS3876.matrix', index_col=0)
        print(df1.shape)
        df1 = quantileNormalize(df1)
        print(df1.shape)
        df1 = percentile(df1, 10)
        print(df1.shape)
        df1 = df_mean_index(df1)
        print(df1.shape)
        df1.to_csv(filePath + fileName + '_qpm.csv')


if __name__ == "__main__":
    sys.exit(main())


