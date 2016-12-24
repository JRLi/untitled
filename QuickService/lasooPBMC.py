import numpy as np
import pandas as pd
from sklearn import datasets, linear_model
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn import metrics
import os

file_path = 'D:/Project/PBMC/logistic_in/'
out_path = 'D:/Project/PBMC/logistic_out/'
file_list = os.listdir(file_path)
print('File_list:', file_list)
# file_path = 'D:/Rworkspace/GSE19301MCLSL'

for file_name in file_list:
    df = pd.read_table(file_path + file_name, index_col=0)
    target = 'Disease'
    cells = df.columns.tolist()
    print(len(cells), cells[0], cells[-1])
    train_score_dict, test_score_dict, test_fscore_dict = {}, {}, {}
    reg = linear_model.Lasso(alpha=0.1)

