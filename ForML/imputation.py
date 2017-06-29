import numpy as np
import pandas as pd
import os
from imputer import Imputer
impute = Imputer()

def openDF(in_path, direct = 'f'):
    fpath, fname = os.path.split(in_path)
    fbase, fext = os.path.splitext(fname)
    df = pd.read_csv(in_path, index_col=0) if fext == '.csv' else pd.read_table(in_path, index_col=0)
    if direct == 't':
        df = df.transpose()
    return df, fbase


def imputation(df_input, col, cat=False, knn_n=10):
    df = df_input.copy()
    imputed_array = impute.knn(X=df, column=col, k=knn_n, is_categorical=cat)
    df_imputed = pd.DataFrame(imputed_array)
    df_imputed.columns = df.columns
    df_imputed.index = df.index
    return df_imputed


def main():
    df1, df_base = openDF('D://Project/RA/RA.txt')
    print(df1.columns)
    print(df1[pd.to_numeric(df1['RF0'], errors='coerce').isnull()])
    df1 = df1.drop('Serial_no2', 1)
    print(df1.shape)

    df1['DAS0'] = df1['DAS0'].replace(-9, np.nan)
    df1['Stag'] = df1['Stag'].replace(-9, np.nan)
    df1['BioDAS1'] = df1['BioDAS1'].replace(-9, np.nan)
    df1['Resp'] = df1['Resp'].replace(9, np.nan)
    df1['RAage'] = df1['RAage'].replace(0, np.nan)
    df1['Bio'] = df1['Bio'].replace(0, np.nan)
    df1['Bioage'] = df1['Bioage'].replace(0, np.nan)
    df_d = df1.dropna()
    df1 = df1.dropna(subset=['Resp', 'Bio'])

    print(df1[df1['Stag'].isnull()])
    print(df1.shape)

    numerical_list = ['RAage', 'Bioage', 'DAS0', 'BioDAS1']
    categorical_list = ['Stag']
    for num_col in numerical_list:
        df1 = imputation(df1, num_col)
    for cat_col in categorical_list:
        df1 = imputation(df1, cat_col, True)

    print(df1.shape)

    decimals = pd.Series([0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 2, 2, 0], index=df1.columns)
    df1 = df1.round(decimals)
    df_d.to_csv('D://Project/RA/RA_drop.txt', sep='\t')
    df1.to_csv('D://Project/RA/RA_imputed.txt', sep='\t')


if __name__ == '__main__':
    main()

