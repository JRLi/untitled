import pandas as pd
from sklearn import datasets, linear_model
from sklearn.model_selection import train_test_split
from sklearn import metrics
import os

file_path = 'D:/Project/PBMC/logistic_in/'
out_path = 'D:/Project/PBMC/logistic_out/'
file_list = os.listdir(file_path)
print('File_list:', file_list)

alpha=0.1

for file_name in file_list:
    print('Processing file:', file_name)
    df = pd.read_table(file_path + file_name, index_col=0)
    target = 'Disease'      # As a 'Y'
    cells = df.columns.tolist()
    features = cells[:-1]   # Modify: Use all features at the same time.
    print(len(features), features[0], features[-1])
    train_score_dict, test_score_dict, test_fscore_dict = {}, {}, {}
    train_df, test_df = train_test_split(df, train_size=0.8, random_state=1)
    lasso_reg = linear_model.Lasso(alpha=alpha)
    lasso_reg.fit(train_df[features], train_df[target])
    test_y_predict = lasso_reg.predict(test_df[features])
    train_y_predict = lasso_reg.predict(train_df[features])
    r2_score_lasso = metrics.r2_score(test_df[target], test_y_predict)
    r2_score_lasso_train = metrics.r2_score(train_df[target], train_y_predict)
    print(test_y_predict)
    print(lasso_reg)
    print('r^2 in test data: %f' % r2_score_lasso)
    print('r^2 in train data: %f' % r2_score_lasso_train)
