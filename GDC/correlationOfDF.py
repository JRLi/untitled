import numpy as np
import pandas as pd
import argparse
use_message = '''
    Need Python3 and numpy, pandas, scipy.
    Usage: correlationOfDF.py df1 df2
'''

def args_parse():
    parser = argparse.ArgumentParser(description=use_message)
    parser.add_argument('profile', nargs='+', help="drug expression profile, separate by space")
    args = parser.parse_args()
    return args

def df_mean_index(df_input):
    # gp1 = df1.groupby(df1.index)
    # df1 = gp1.mean()
    gp = df_input.groupby(df_input.index)
    return gp.mean()

df9 = pd.read_table('D://Project/drs/forTest9x7.txt', index_col= 0)
df8 = pd.read_table('D://Project/drs/forTest8x8.txt', index_col= 0)     # drug
df8 = df_mean_index(df8)    # mean of drug
print(df8)
print('aaaaa')
print(df8.corrwith(df9.iloc[:, 1], axis=0))
print('bbb')
# print(df12.corrwith(df9.iloc[:, 1], axis=1))
# print(df8.corrwith(df9.iloc[:, 0]))
# print(df9.columns[0])
df2 = pd.DataFrame(columns=df9.columns, index=df8.columns)
df5c = pd.DataFrame(columns=df9.columns, index=df8.columns)
df5p = pd.DataFrame(columns=df9.columns, index=df8.columns)
df3 = df9.merge(df8,left_index=True, right_index=True, how='left').dropna()
print('merge test')
print(df3)

# df4 = df3.corr(method='pearson', min_periods=1)
for i in range(len(df9.columns)):
    df2.iloc[:, i] = df8.corrwith(df9.iloc[:, i])
print(len(df3.columns))
print(df2)

for i in range(len(df9.columns)):
    c1 = df3.columns[i]
    c_list, p_list = [], []
    for j in range(len(df9.columns), len(df3.columns)):
        c2 = df3.columns[j]
        ftols = pd.ols(y = df3[c1], x = df3[c2], intercept=True)
        c_list.append(df3[c1].corr(df3[c2]))
        p_list.append(ftols.f_stat['p-value'])
    df5c[c1] = pd.Series(c_list).values
    df5p[c1] = np.array(p_list)
print(df5c)
print('')
print(df5p)
