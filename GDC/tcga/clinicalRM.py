import numpy as np
import pandas as pd

clinical_path = 'D:/Project/drs/GDC/all_clin_XML.csv'

df1 = pd.read_csv(clinical_path)
print(df1.shape)
df1 = df1.drop_duplicates()
print(df1.shape)
df1 = df1.dropna(axis = 1, how='all')
print(df1.shape)

df1.to_csv('D:/Project/drs/GDC/all_clin_XML_processed.csv', index=False)