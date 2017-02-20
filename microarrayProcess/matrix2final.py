from microarrayProcess.arrayDFquantilePercentile import df_mean_columns, hgnc_dict
import pandas as pd
filePath = 'E:/StringTemp/GDS/'
fileName = 'GDS3876.matrix_qpm.csv'
hgncFile = 'hgnc_complete_set.txt'
df = pd.read_csv(filePath + fileName, index_col= 0)
print(df.shape)
# print(df.columns)
df.columns = ['LNT2D', 'LNT2D', 'LNT2D', 'LNT2D', 'LNT2D', 'ONT2D', 'ONT2D', 'ONT2D', 'ONT2D', 'OWT2D', 'OWT2D',
              'OWT2D', 'OWT2D', 'OWT2D', 'OPT2D', 'OPT2D', 'OPT2D', 'OPT2D']
print(df.shape, df.columns, sep='\n')
df = df_mean_columns(df)
print(df.iloc[0:5,])

symbol2HgncDict = hgnc_dict(filePath + hgncFile)
print(len(symbol2HgncDict))
