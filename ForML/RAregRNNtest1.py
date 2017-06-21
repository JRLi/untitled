from keras.models import model_from_json, Sequential
from keras.utils.np_utils import to_categorical
import numpy as np
import pandas as pd
import os, sys
import scipy.stats as stats
import tensorflow as tf
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, f1_score, roc_curve, auc
import matplotlib.pyplot as plt
from keras.layers import Dense, Dropout
from keras.wrappers.scikit_learn import KerasRegressor
from tensorflow.python.client import device_lib
import keras

print(device_lib.list_local_devices())
print("python:{}, keras:{}, tensorflow: {}".format(sys.version, keras.__version__, tf.__version__))

in_path = 'D:/Project/RA/in/'
out_path = 'D:/Project/RA/out/'
file = 'RA_drop.txt'
fpr, tpr, roc_auc = {}, {}, {}
output_dim1, output_dim2, batch_size, nb_epoch = 8, 5, 10, 100


def parse_label(label_type_list, label_type):
    types = list(set(label_type_list))
    types.sort()
    types = {types[i]: i for i in range(len(types))}
    if label_type == 'c':
        labels = [[types[label_type]] for label_type in label_type_list]
        labels = np.array(labels)
        labels = to_categorical(labels, len(types))
    else:
        labels = label_type_list
    print(labels)
    return labels, types


def get_data(path, phenotype, type='c', proportion=0.8, mode='n'):
    fpath, fname = os.path.split(path)
    fbase, fext = os.path.splitext(fname)
    df = pd.read_csv(path, index_col=0) if fext == '.csv' else pd.read_table(path, index_col=0)
    target = phenotype  # As a 'Y', target
    if mode != 'n':   # mode 1 indicate feature list is first row, and sample is first column
        df = df.transpose

    train_df, test_df = train_test_split(df, train_size=proportion, random_state=1)  # train_size
    cells = df.columns.tolist()
    cells.remove(target)
    X_train, y_train = train_df[cells].values.astype('float32'), train_df[target]
    X_test, y_test = test_df[cells].values.astype('float32'), test_df[target]
    # np.linalg.norm
    labels_train, types_train = parse_label(y_train.tolist(), type)
    labels_test, types_test = parse_label(y_test.tolist(), type)
    return X_train, labels_train, types_train, X_test, labels_test, types_test


def corr(a, b):
    return stats.pearsonr(a, b)


def score(true_y, predict_y, mode=1):
    accuracy = accuracy_score(true_y, predict_y)
    if mode == 1:
        f1score = f1_score(true_y, predict_y)
    elif mode == 2:
        f1score = f1_score(true_y, predict_y, average='macro')
    else:
        f1score = 'No calculate F1 score.'  # mode == 1 indicate binary data
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


def roc_plot(fpr_dict, tpr_dict, roc_dict, out_path, title = 'keras'):
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
    plt.title('Receiver operating characteristic of ' +title)
    plt.legend(loc="lower right")
    fig1 = plt.gcf()
    # fig1.set_size_inches(6, 6)
    # fig1.set_dpi(300)
    fig1.savefig(out_path + title + '.png')
    plt.show()


def main():
    # if the data format is not standard, must assign the mode parameter to not 1.
    data_train, labels_train, types_train, data_test, labels_test, types_test = get_data(in_path + file, 'Resp', 'n', 0.67)
    print('type(data_train):', type(data_train))
    print(data_train.shape)
    print(types_train)
    # print('data_train.iloc[1]:',data_train.iloc[2])

    np.random.seed(1)
    model = Sequential()
    model.add(Dense(units=output_dim1, input_dim=len(data_train[0]), activation="relu"))
    model.add(Dropout(0.2))
    model.add(Dense(units=output_dim2, activation="relu"))
    model.add(Dropout(0.3))
    # model.add(Dense(len(types_train), activation="sigmoid"))
    model.add(Dense(1, activation="sigmoid"))
    model.compile(loss="mse", optimizer="rmsprop")
    print(model.summary())

    model.fit(data_train, labels_train, epochs=nb_epoch, batch_size=batch_size, verbose=1, validation_split=0.1)
    train_score = model.predict(data_train)
    predict_types = model.predict_classes(data_train)

    print(train_score.shape, predict_types.shape)
    print(train_score[0:3, :])
    print(predict_types[0:3])

    fpr['Training'], tpr['Training'], roc_auc['Training'] = generate_results(labels_train[:, 0], train_score[:, 0])
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

    roc_plot(fpr, tpr, roc_auc, out_path, file)

if __name__ == '__main__':
    main()