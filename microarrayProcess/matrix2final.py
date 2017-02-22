from microarrayProcess.arrayDFquantilePercentile import df_mean_columns, hgnc_dict
import pandas as pd

filePath = 'E:/StringTemp/GDS/'
fileName = 'GDS3876.matrix_qpm.csv'
hgncFile = 'hgnc_complete_set.txt'

df = pd.read_csv(filePath + fileName, index_col= 0)
print(df.shape)

df.columns = ['LNT2D', 'LNT2D', 'LNT2D', 'LNT2D', 'LNT2D', 'ONT2D', 'ONT2D', 'ONT2D', 'ONT2D', 'OWT2D', 'OWT2D',
              'OWT2D', 'OWT2D', 'OWT2D', 'OPT2D', 'OPT2D', 'OPT2D', 'OPT2D']
print(df.shape, df.columns, sep='\n')

df = df_mean_columns(df)    # fixed step

# debug step
print(df.iloc[0:5,])
print(df.iloc[0:5,0]/df.iloc[0:5,1])
print(df['LNT2D'].iloc[0:5,] / df['ONT2D'].iloc[0:5,])

symbol_list = df.index.tolist()      # fixed step
symbol2HgncDict = hgnc_dict(filePath + hgncFile)    # fixed step
print(len(symbol2HgncDict), len(set(symbol2HgncDict.keys())), len(set(symbol2HgncDict.values())))

hgnc_list = [symbol2HgncDict.get(symbol,'-') for symbol in symbol_list]      # fixed step

out_file = ['ONT2DvsLNT2D', 'OWT2DvsONT2D', 'OPT2DvsOWT2D', 'OPT2DvsONT2D', 'OPT2DvsLNT2D', 'OWT2DvsLNT2D']

for versus in out_file:
    case, control = versus.split('vs')
    df2 = pd.DataFrame(columns=['case', 'control', 'hgnc', 'sym', 'fc'])
    df2['sym'] = pd.Series(symbol_list).values
    df2['hgnc'] = pd.Series(hgnc_list).values
    df2['case'] = case
    df2['control'] = control
    df2['fc'] = (df[case] / df[control]).values
    df2.to_csv(filePath + versus + '.csv', header=False, index=False)
