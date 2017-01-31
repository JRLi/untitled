from keras.models import model_from_json
from keras.utils.np_utils import to_categorical
import numpy as np
import pandas as pd
import scipy.stats as stats
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, f1_score, roc_curve, auc
import matplotlib.pyplot as plt


def parse_label(label_type_list):
    types = list(set(label_type_list))
    types.sort()
    types = {types[i]: i for i in range(len(types))}
    labels = [[types[label_type]] for label_type in label_type_list]
    labels = np.array(labels)
    labels = to_categorical(labels, len(types))
    return labels, types


def get_data(path, phenotype, proportion, mode=1):
    df = pd.read_table(path, index_col=0)
    target = phenotype  # As a 'Y', target
    if mode != 1:   # mode 1 indicate feature list is first row, and sample is first column
        df = df.transpose
    train_df, test_df = train_test_split(df, train_size=proportion, random_state=1)  # train_size
    cells = df.columns.tolist()
    # features = cells[:-1]
    cells.remove(target)
    data_train, target_train = train_df[cells].values.astype('float32'), train_df[target]
    data_test, target_test = test_df[cells].values.astype('float32'), test_df[target]
    # np.linalg.norm
    labels_train, types_train = parse_label(target_train.tolist())
    labels_test, types_test = parse_label(target_test.tolist())
    return data_train, labels_train, types_train, data_test, labels_test, types_test


def corr(a, b):
    return stats.pearsonr(a, b)


def score(true_y, predict_y, mode=1):
    accuracy = accuracy_score(true_y, predict_y)
    f1score = f1_score(true_y, predict_y) if mode == 1 else 'No calculate F1 score.'  # mode == 1 indicate binary data
    return accuracy, f1score


def save_model(model, path, file_name):
    json_string = model.to_json()
    with open(path + file_name + ".model.json", "w") as output:
        output.writelines(json_string)
    model.save_weights(path + file_name+ ".weights.h5")


def load_model(path, file_name):
    with open(path + file_name + ".model.json", "r") as input:
        model = model_from_json(input.readline())
    model.load_weights(path + file_name + ".weights.h5")
    return model


def generate_results(y_test, y_score):
    fpr, tpr, _ = roc_curve(y_test, y_score)
    roc_auc = auc(fpr, tpr)
    return fpr, tpr, roc_auc


def roc_plot(fpr_dict, tpr_dict, roc_dict):
    plt.figure(1)
    colors = ['darkgreen', 'darkred']
    lw = 2
    for i, color in zip(fpr_dict.keys(), colors):
        plt.plot(fpr_dict[i], tpr_dict[i], color=color, lw=lw,
                 label='ROC curve of {0} (area = {1:0.2f})'
                       ''.format(i, roc_dict[i]))
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic of keras')
    plt.legend(loc="lower right")
    plt.show()
