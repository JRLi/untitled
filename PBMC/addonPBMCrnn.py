from keras.models import model_from_json
from keras.utils.np_utils import to_categorical
from numpy import transpose
import numpy as np
import pandas as pd
import scipy.stats as stats
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, f1_score


def parse_label(label_type_list):
    types = list(set(label_type_list))
    types.sort()
    types = {types[i]: i for i in range(len(types))}
    labels = [[types[lable_type]] for lable_type in label_type_list]
    labels = np.array(labels)
    labels = to_categorical(labels, len(types))
    return labels, types


def get_data(path, phenotype, proportion, mode=1):
    df = pd.read_table(path, index_col=0)
    target = phenotype  # As a 'Y', target
    if mode != 1:
        df = df.transpose
    cells = df.columns.tolist()
    features = cells[:-1]
    train_df, test_df = train_test_split(df, train_size=proportion, random_state=1)    # train_size
    data_train, target_train = train_df[features].values.astype('float32'), train_df[target]
    data_test, target_test = test_df[features].values.astype('float32'), test_df[target]
    # np.linalg.norm
    labels_train, types_train = parse_label(target_train.tolist())
    labels_test, types_test = parse_label(target_test.tolist())
    return data_train, labels_train, types_train, data_test, labels_test, types_test


def corr(a, b):
    return stats.pearsonr(a, b)[0]


def score(true_y, predict_y, mode = 1):
    accuracy = accuracy_score(true_y, predict_y)
    f1score = f1_score(true_y, predict_y) if mode == 1 else 'No calculate F1 score.'
    return accuracy, f1score


def save_model(model, path):
    json_string = model.to_json()
    with open(path + "model.json", "w") as output:
        output.writelines(json_string)
    model.save_weights(path + "weights.h5")


def load_model(path):
    with open(path + "model.json", "r") as input:
        model = model_from_json(input.readline())
    model.load_weights(path + "weights.h5")
    return model
