#!/usr/bin/env python3.5
"""
Created by JRLi on 2016/12/06 for python learning
"""
import sys
import argparse
import os
import numpy as np
import pandas as pd
import scipy.stats as stats
from keras.models import Sequential, model_from_json
from keras.layers import Dense, Dropout
import tensorflow as tf
from keras.utils.np_utils import to_categorical
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, f1_score, roc_curve, auc
import matplotlib.pyplot as plt

tf.python.control_flow_ops = tf

use_message = '''
    Usage:  fastaqcr [-r] fa [-b BASE] [-i] input [input ...]
            fastaqcr [-r] fq [-p Left Right] [-s [SINGLE]]
    For help, please use "fastaqcr fa -h" or "fastaqcr fq -h".
    To calculate NGS reads length and the number of reads.
    To perform the reverse complement of fasta/fastq file.
    Uniform fasta file line length: one line or 60mer/per line.
'''

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


def args_parse():
    pass


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
    colors = ['darkgreen', 'darkred', 'darkblue', 'aqua', 'darkorange', 'pink', 'purple', 'black']  # max colors are 8
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


def main(argv=None):
    try:
        if argv is None:
            argv = args_parse()
            in_path = 'D:/Project/PBMC/logistic_in/'
            out_path = 'D:/Project/PBMC/logistic_out/'
            file = 'GSE16129_19301_0_SAI_1_exacerbationAsthma'
            fpr = dict()
            tpr = dict()
            roc_auc = dict()

            output_dim1 = 64
            output_dim2 = 32
            batch_size = 10
            nb_epoch = 1000

            # if the data format is not standard, must assign the mode parameter to not 1.
            data_train, labels_train, types_train, data_test, labels_test, types_test = get_data(in_path + file,
                                                                                                 'Disease', 0.8)
            print('type(data_train):', type(data_train))
            print(data_train.shape)
            # print('data_train.iloc[1]:',data_train.iloc[2])

            model = Sequential()

            model.add(Dense(output_dim=output_dim1, input_dim=len(data_train[0]), activation="relu"))
            model.add(Dropout(0.25))
            model.add(Dense(output_dim2, activation="relu"))
            model.add(Dropout(0.5))
            model.add(Dense(len(types_train), activation="sigmoid"))

            model.compile(loss="mse", optimizer="rmsprop")

            model.fit(data_train, labels_train, nb_epoch=nb_epoch, batch_size=batch_size)

            train_score = model.predict(data_train)
            predict_types = model.predict_classes(data_train)

            fpr['Training'], tpr['Training'], roc_auc['Training'] = generate_results(labels_train[:, 0],
                                                                                     train_score[:, 0])
            correlation, p_value = corr(labels_train.argmax(1).tolist(), predict_types.tolist())
            # if the labels are not binary type, must assign the mode parameter to not 1.
            accuracy, f1score = score(labels_train.argmax(1).tolist(), predict_types.tolist())
            print("\nFor Training\nTypes:", types_train)
            print("True Types:", labels_train.argmax(1))
            print("Predict Types:", predict_types)
            print("Corr: {}\np-value: {}".format(correlation, p_value))
            print('Accuracy: {}\nF1 score: {}'.format(accuracy, f1score))

            save_model(model, out_path, file)
            model = load_model(out_path, file)

            test_score = model.predict(data_test)
            predict_types = model.predict_classes(data_test)
            fpr['Testing'], tpr['Testing'], roc_auc['Testing'] = generate_results(labels_test[:, 0], test_score[:, 0])
            correlation, p_value = corr(labels_test.argmax(1).tolist(), predict_types.tolist())
            accuracy, f1score = score(labels_test.argmax(1).tolist(), predict_types.tolist())
            print("\nFor Testing\nTypes:", types_test)
            print("True Types:", labels_test.argmax(1))
            print("Predict Types:", predict_types)
            print("Corr: {}\np-value: {}".format(correlation, p_value))
            print('Accuracy: {}\nF1 score: {}'.format(accuracy, f1score))

            roc_plot(fpr, tpr, roc_auc)
        print('Done.')

    except Usage as err:
        print(sys.stderr, err.msg)
        print(sys.stderr, "Terminated, for help use -h or --help")
        return 2

if __name__ == "__main__":
    sys.exit(main())
