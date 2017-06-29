#!/usr/bin/env python3
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn import preprocessing
from sklearn.model_selection import train_test_split
from keras.models import Sequential
from keras.layers import Dense, Dropout

np.random.seed(10)  # set random seed
in_path = 'D:/Project/RA/in/'
out_path = 'D:/Project/RA/out/'
file = 'RA.txt'


def openDF(in_path, direct='n'):
    fpath, fname = os.path.split(in_path)
    fbase, fext = os.path.splitext(fname)
    df = pd.read_csv(in_path, index_col=0) if fext == '.csv' else pd.read_table(in_path, index_col=0)
    if direct == 't':
        df = df.transpose()
    return df


def preprocess_data(raw_df):
    df = raw_df.copy()
    df = num2nan(df)
    df = df.dropna(subset=['Resp', 'Bio'])
    df['drugTime'] = df['Bioage'] - df['RAage']
    df['drugTime'][df['drugTime'] < 0] = 0.0
    df = df.drop(['Bioage', 'RAage'], axis=1)
    all_mean = df.mean()
    print(df.isnull().sum())
    df['DAS0'] = df['DAS0'].fillna(all_mean.DAS0)
    df['Stag'] = df['Stag'].fillna(all_mean.Stag)
    df['drugTime'] = df['drugTime'].fillna(all_mean['drugTime'])
    df['Resp'] = df['Resp'].map({2.0: 1.0, 1.0: 0.0, 0.0: 0.0})
    #df['Resp'] = df['Resp'].replace(2.0, 1.0)
    cols_list = df.columns.tolist()  # must need tolist()
    cols_list = cols_list[-1:] + cols_list[:-1]
    df = df[cols_list]
    label = df.values[:, -1]
    features = df.values[:, 0:-1]
    print(df[:2])
    print(features[:2])
    print(label[:2])
    print('label_shape:', label.shape, 'features_shape:', features.shape)
    # standardize scale for features
    minmax_scale = preprocessing.MinMaxScaler(feature_range=(0, 1))
    scaledFeathres = minmax_scale.fit_transform(features)
    # return feature and label
    return scaledFeathres, label


def num2nan(df_in):
    df_in['DAS0'] = df_in['DAS0'].replace(-9, np.nan)
    df_in['Stag'] = df_in['Stag'].replace(-9, np.nan)
    df_in['Resp'] = df_in['Resp'].replace(9, np.nan)
    df_in['RAage'] = df_in['RAage'].replace(0, np.nan)
    df_in['Bio'] = df_in['Bio'].replace(0, np.nan)
    df_in['Bioage'] = df_in['Bioage'].replace(0, np.nan)
    return df_in


def mlp(X, y):
    model = Sequential()
    model.add(Dense(units=40, input_dim=X.shape[1], kernel_initializer='uniform', activation='relu'))
    model.add(Dropout(0.2))
    model.add(Dense(units=30, kernel_initializer='uniform', activation='relu'))
    model.add(Dropout(0.1))
    model.add(Dense(units=1, kernel_initializer='uniform', activation='sigmoid'))
    print(model.summary())
    model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
    train_history = model.fit(x=X, y=y, validation_split=0.3, epochs=30, batch_size=10, verbose=2)
    return train_history, model


def show_train_history(train_history, train, validation):
    plt.plot(train_history.history[train])
    plt.plot(train_history.history[validation])
    plt.title('Train History')
    plt.ylabel(train)
    plt.xlabel('Epoch')
    plt.legend(['train', 'validation'], loc='upper left')
    plt.show()


def main():
    all_df = openDF(in_path + file)
    cols = all_df.columns[[1, 2, 3, 4, 5, 7, 9, 10, 11, 13]]
    print(len(cols), cols)
    all_df = all_df[cols]
    train_df, test_df = train_test_split(all_df, train_size=0.7, random_state=1)  # train_size
    train_features, train_label = preprocess_data(train_df)
    test_features, test_label = preprocess_data(test_df)
    print(train_features[:2])
    print(train_label[:2])
    print(test_features[:2])
    print(test_label[:2])
    # training model
    train_history, model = mlp(train_features, train_label)
    show_train_history(train_history, 'acc', 'val_acc')
    show_train_history(train_history, 'loss', 'val_loss')
    scores = model.evaluate(x=test_features, y=test_label)
    print('accuracy:', scores[1])
    model.save_weights(out_path + 'mlp_ra.h5')

    all_features, all_label = preprocess_data(all_df)
    all_probability = model.predict(all_features)
    all_class_p = model.predict_classes(all_features)
    #print(all_probability[all_probability < 0.5])
    print(all_probability[:10])
    print(all_label[:10])
    print(all_class_p[:10])

if __name__ == '__main__':
    main()