import pandas as pd
from sklearn import datasets, linear_model
import matplotlib.pyplot as plt
from sklearn import metrics
from sklearn.model_selection import train_test_split
import os

file_path = 'D:/Project/PBMC/logistic_in/'
out_path = 'D:/Project/PBMC/logistic_out/'
file_list = os.listdir(file_path)
print(file_list)
# file_path = 'D:/Rworkspace/GSE19301MCLSL'

for file_name in file_list:
    df = pd.read_table(file_path + file_name, index_col=0)
    target = 'Disease'
    cells = df.columns.tolist()
    train_score_dict = {}
    test_score_dict = {}
    # for i in range(0, 5):
    features = ['T.8NVE.SP.OT1', 'GN.ARTH.BM', 'NKT.4pos.LV', 'GN.BM', 'TGD.VG2pos.SP']   # can use for loop to do one cell iteration
    print(''.join(features))
    train_df, test_df = train_test_split(df, train_size=0.8, random_state=1)
    # print(train_df.shape, test_df.shape)
    # print(train_df[features].shape, train_df[target].shape)
    regularization = 1.0    # 1e5
    logreg = linear_model.LogisticRegression(C=regularization)
    logreg.fit(train_df[features], train_df[target])
    preds = logreg.predict_proba(test_df[features])
    preds_1 = preds[:,1]
    fpr, tpr, _ = metrics.roc_curve(test_df[target], preds_1)
    auc = metrics.auc(fpr, tpr)

    print(logreg.predict(test_df[features]).tolist(), test_df[target].tolist(), sep="\n")
    score_Train = logreg.score(train_df[features], train_df[target])
    score_Test = logreg.score(test_df[features], test_df[target])
    print(features, score_Train, score_Test)
    train_score_dict[''.join(features)] = score_Train
    test_score_dict[''.join(features)] = score_Test

    plt.figure()
    lw = 2
    plt.plot(fpr, tpr, color='darkorange', lw = lw, label='ROC curve (area = %0.2f)' % auc)
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic of 5 cells logistic')
    plt.legend(loc="lower right")
    plt.show()
