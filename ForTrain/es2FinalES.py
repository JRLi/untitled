#!/usr/bin/env python
import pandas as pd
import sys
import os
import datetime

def final_ES(df_input):
    df_up = df_input.iloc[:, 0:len(df_input.columns) // 2]
    df_dn = df_input.iloc[:, len(df_input.columns) // 2:len(df_input.columns)]
    print(df_up.shape)
    print(df_dn.shape)
    df_all = df_up - df_dn.values
    print(df_all.shape)
    df_all.columns = [c[:-6] for c in df_all.columns]
    return df_all

st = datetime.datetime.now()
df1 = pd.read_csv(sys.argv[1], index_col=0)
fpath, fname = os.path.split(sys.argv[1])
fbase, fext = os.path.splitext(fname)
df2 = final_ES(df1)
df2.to_csv('./{}_finalES{}'.format(fbase, fext), na_rep='NA')
td = datetime.datetime.now() - st
print("\t[Info] Spending time={0}!".format(td))