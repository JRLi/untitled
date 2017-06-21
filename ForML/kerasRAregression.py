from keras.models import Sequential
from keras.layers.core import Dense, Activation
import pandas as pd
import numpy as np
import os, sys
from tensorflow.python.client import device_lib
from keras.utils.np_utils import to_categorical
from sklearn.model_selection import train_test_split
import keras
import tensorflow as tf
print(device_lib.list_local_devices())

print("python:{}, keras:{}, tensorflow: {}".format(sys.version, keras.__version__, tf.__version__))
in_path = 'D:/Project/RA/in/'
out_path = 'D:/Project/RA/out/'
file = 'RA_drop.txt'


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


X_train, y_train, types_train, X_test, y_test, types_test = get_data(in_path + file, 'Resp', 'n', 0.67)
print(len(X_train[0]))

model = Sequential()
model.add(Dense(2, input_dim=len(X_train[0]),init='uniform', activation='linear'))
model.compile(loss='mse', optimizer='rmsprop')

model.fit(X_train, y_train, epochs=100, batch_size=16,verbose=0)
model.fit(X_train, y_train, epochs=1, batch_size=16,verbose=1)
score = model.evaluate(X_test, y_test, batch_size=16)