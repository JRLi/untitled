import numpy as np
import pandas as pd
from sklearn import datasets, linear_model
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn import metrics
from itertools import cycle
import os

file_path = 'D:/Project/PBMC/logistic_in/'
out_path = 'D:/Project/PBMC/logistic_out/'
file_list = os.listdir(file_path)
print(file_list)
# file_path = 'D:/Rworkspace/GSE19301MCLSL'
fpr = dict()
tpr = dict()
roc_auc = dict()

for file_name in file_list:
    df = pd.read_table(file_path + file_name, index_col=0)
    target = 'Disease'
    cells = df.columns.tolist()
    train_score_dict = {}
    test_score_dict, test_fscore_dict = {}, {}
    for i in range(0, len(cells) - 1):
        # for i in range(0, 5):
        features = cells[i:i+1]   # can use for loop to do one cell iteration
        print(''.join(features))
        train_df, test_df = train_test_split(df, train_size=0.8, random_state=1)
        # print(train_df.shape, test_df.shape)
        regularization = 1.0    # 1e5
        logreg = linear_model.LogisticRegression(C=regularization)
        logreg.fit(train_df[features], train_df[target])
        preds = logreg.predict_proba(test_df[features])
        preds_1 = preds[:, 1]
        cell = ''.join(features)
        fpr[cell], tpr[cell], _ = metrics.roc_curve(test_df[target], preds_1)
        roc_auc[cell] = metrics.auc(fpr[cell], tpr[cell])

        p_list = logreg.predict(test_df[features]).tolist()
        t_list = test_df[target].tolist()
        print(p_list, t_list, sep="\n")
        score_Train = logreg.score(train_df[features], train_df[target])
        score_Test = logreg.score(test_df[features], test_df[target])
        # scoreTestA = metrics.accuracy_score(t_list,p_list)   # fraction of correctly classified samples
        scoreTestF = metrics.f1_score(t_list, p_list)
        print(features, score_Train, score_Test, scoreTestF)    # 3rd score is metrics score
        # train_score_dict[''.join(features)] = score_Train
        test_score_dict[''.join(features)] = score_Test
        test_fscore_dict[''.join(features)] = scoreTestF
        # print(test_df[features][:2])
        # print(logreg.predict_proba(test_df[features][:2]))
    test_series = pd.Series(test_score_dict).sort_values()
    test_fseries = pd.Series(test_fscore_dict).sort_values()
    # print(test_series)
    print(len(cells), cells[-1])
    print(len(test_series), test_series.iloc[-5:], sep="\n")
    print(len(test_fseries), test_fseries.iloc[-5:], sep="\n")
    ser_a = test_fseries.iloc[-5:]
    for i in ser_a.index:
        print(i)
    plt.figure(1)
    plt.plot(test_series.values)
    plt.xlabel('Sorted immune cells')
    plt.ylabel('Mean accuracy')
    plt.title(file_name + '_accuracy')
    plt.figure(2)
    plt.plot(test_fseries.values)
    plt.xlabel('Sorted immune cells')
    plt.ylabel('Mean accuracy')
    plt.title(file_name + '_f-score')
    plt.figure(3)
    colors = cycle(['aqua', 'darkorange', 'cornflowerblue', 'green', 'red'])
    lw = 2
    for i, color in zip(ser_a.index, colors):
        plt.plot(fpr[i], tpr[i], color=color, lw=lw,
             label='ROC curve of class {0} (area = {1:0.2f})'
             ''.format(i, roc_auc[i]))
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic of logistic regression: ' + file_name)
    plt.legend(loc="lower right")
    plt.show()
