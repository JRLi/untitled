#!/usr/bin/env python3
import urllib.request
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn import preprocessing
from keras.models import Sequential
from keras.layers import Dense, Dropout

np.random.seed(10)  # set random seed
url = 'http://biostat.mc.vanderbilt.edu/wiki/pub/Main/DataSets/titanic3.xls'
filepath = 'E://StringTemp/titanic3.xls'


def preprocess_data(raw_df):
    # 1. remove name
    #df = raw_df.copy()
    df =raw_df.drop(['name'], axis=1)  # drop name, drop can input string or list, axis = 1 indicate row-wise
    # 2. find nan columns
    print(df.isnull().sum())  # age, fare, embarked have nan
    # 3. get mean and fill nan with mean
    af_mean = df[['age', 'fare']].mean()  # get mean of age and fare
    df['age'] = df['age'].fillna(af_mean.age)  # fill nan with age mean
    df['fare'] = df['fare'].fillna(af_mean['fare'])  # fill nan with fare mean
    # 4. transfer sex to binary
    df['sex'] = df['sex'].map({'female': 0, 'male': 1}).astype(int)
    # 5. coding embarked to OneHot encoding
    x_OneHot_df = pd.get_dummies(df, columns=['embarked'])
    print('onehot_df_shape:', x_OneHot_df.shape)
    # training needs ndarray format
    label = x_OneHot_df.values[:, 0]
    features = x_OneHot_df.values[:, 1:]
    print('label_shape:', label.shape, 'features_shape:', features.shape)
    # standardize scale for features
    minmax_scale = preprocessing.MinMaxScaler(feature_range=(0, 1))
    scaledFeathres = minmax_scale.fit_transform(features)
    # return feature and label
    return scaledFeathres, label


def show_train_history(train_history, train, validation):
    plt.plot(train_history.history[train])
    plt.plot(train_history.history[validation])
    plt.title('Train History')
    plt.ylabel(train)
    plt.xlabel('Epoch')
    plt.legend(['train', 'validation'], loc='upper left')
    plt.show()


def mlp(X, y):
    model = Sequential()
    model.add(Dense(units=40, input_dim=X.shape[1], kernel_initializer='uniform', activation='relu'))
    model.add(Dropout(0.2))
    model.add(Dense(units=30, kernel_initializer='uniform', activation='relu'))
    model.add(Dropout(0.1))
    model.add(Dense(units=1, kernel_initializer='uniform', activation='sigmoid'))
    print(model.summary())
    model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
    train_history = model.fit(x=X, y=y, validation_split=0.1, epochs=30, batch_size=30, verbose=2)
    return train_history, model


def main():
    # if data not exist, download data
    if not os.path.isfile(filepath):
        result = urllib.request.urlretrieve(url, filepath)
        print('downloaded:', result)
    # get useful information
    all_df = pd.read_excel(filepath)    # read data frame from excel
    cols = all_df.columns[[1, 2, 0, 3, 4, 5, 6, 8, 10]]     # get useful columns
    print(cols)
    all_df = all_df[cols]
    # get random training and testing
    msk = np.random.rand(len(all_df)) < 0.8     # random 0 ~ 1 values of the given shape
    train_df = all_df[msk]  # msk is a boolean array with len(all_df) shape
    test_df = all_df[~msk]  # reverse boolean with msk
    print('all, train, test:', len(all_df), len(train_df), len(test_df))
    # data pre-process
    train_features, train_label = preprocess_data(train_df)
    test_features, test_label = preprocess_data(test_df)
    print(train_features.shape[1])
    # training model
    train_history, model = mlp(train_features, train_label)
    show_train_history(train_history, 'acc', 'val_acc')
    show_train_history(train_history, 'loss', 'val_loss')
    scores = model.evaluate(x=test_features, y=test_label)
    print('accuracy:', scores[1])
    model.save_weights('E://StringTemp/mlp_titanic.h5')

    # build new df with jack and rose
    Jack = pd.Series([0, 'Jack', 3, 'male', 23, 1, 0, 5.000, 'S'])
    Rose = pd.Series([1, 'Rose', 1, 'female', 20, 1, 0, 100.0000, 'S'])
    JR_df = pd.DataFrame([list(Jack), list(Rose)], columns=all_df.columns)
    all_df = pd.concat([all_df, JR_df])
    print(all_df[-2:])

    all_features, all_labels = preprocess_data(all_df)
    all_probability = model.predict(all_features)
    print(all_probability[:10])

    # integrate all_df and all_probability
    all_df_withP = all_df.copy()
    all_df_withP.insert(len(all_df.columns), 'probability', all_probability)
    print(all_df_withP[-2:])
    df_p90n1 = all_df_withP[(all_df_withP['survived'] == 0) & (all_df_withP['probability'] > 0.9)]
    print(df_p90n1)

if __name__ == '__main__':
    main()
