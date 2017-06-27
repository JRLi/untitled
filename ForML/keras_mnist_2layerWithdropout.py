#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from keras.utils import np_utils
from keras.datasets import mnist
from keras.models import Sequential
from keras.layers import Dense
from keras.layers import Dropout

np.random.seed(10)


def plot_image(image):
    fig = plt.gcf()
    fig.set_size_inches(2, 2)
    plt.imshow(image, cmap='binary')
    plt.show()


def plot_image_label_prediction(images, labels, prediction, idx, num=10):
    fig=plt.gcf()
    fig.set_size_inches(12, 14)
    num = 25 if num > 25 else num   # set sub-graph number limitation
    for i in range(0, num):
        ax = plt.subplot(5, 5, 1 + i)   # build 5 x 5 sub-graph
        ax.imshow(images[idx], cmap='binary')   # plot sub-graph of idx
        title = 'label=' + str(labels[idx])     # set title of each sub-graph
        if len(prediction) > 0:     # prediction can be an empty list, indict no prediction
            title += ',predict=' + str(prediction[idx]) # if import prediction result, set title with prediction result
        ax.set_title(title, fontsize=10)    # set title size
        ax.set_xticks([])   # remove the tick of x axis of sub-graph
        ax.set_yticks([])   # remove the tick of y axis of sub-graph
        idx += 1    # next index
    plt.show()


def mlp(X, y):
    model = Sequential()
    model.add(Dense(units=1000, input_dim=len(X[0]), kernel_initializer='normal', activation="relu")) # hidden layer1
    model.add(Dropout(0.5))
    model.add(Dense(units=1000, kernel_initializer='normal', activation="relu"))    # hidden layer 2, 1000
    model.add(Dropout(0.5))
    model.add(Dense(units=y.shape[1], kernel_initializer='normal', activation='softmax'))
    print(model.summary())
    model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
    train_history = model.fit(x=X, y=y, validation_split=0.2, epochs=10, batch_size=200, verbose=2)
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
    # load data
    (X_train_image, y_train_label), (X_test_image, y_test_label) = mnist.load_data()
    print(type(X_test_image))
    print('X_train_image:', X_train_image.shape)
    print('y_train_label:',y_train_label.shape)
    # test plot 1
    #plot_image(X_train_image[1])
    # test plot 2
    #plot_image_label_prediction(X_train_image, y_train_label, [], 0, 10)
    # test plot 3
    #plot_image_label_prediction(X_test_image, y_test_label, [], 0, 10)

    # reshape 28 x 28 array to a 768 x 1 vector
    X_train = X_train_image.reshape(60000, 784).astype('float32')
    X_test = X_test_image.reshape(10000, 784).astype('float32')
    X_train_normalize = X_train / 255
    X_test_normalize = X_test / 255

    # One-hot encoding for labels
    print(y_train_label[:5])
    y_train_OneHot = np_utils.to_categorical(y_train_label)
    y_test_OneHot = np_utils.to_categorical(y_test_label)
    print(y_train_OneHot[:5])
    print(y_train_OneHot.shape)

    train_history, model = mlp(X_train_normalize, y_train_OneHot)
    print(train_history.history)
    show_train_history(train_history, 'acc', 'val_acc')
    show_train_history(train_history, 'loss', 'val_loss')

    scores = model.evaluate(X_test_normalize, y_test_OneHot)
    print(scores)
    print('accuracy=', scores[1])

    prediction = model.predict_classes(X_test)
    print(prediction)
    plot_image_label_prediction(X_test_image, y_test_label, prediction, idx=340, num=25)

    # confusion matrix construction
    cx = pd.crosstab(y_test_label, prediction, rownames=['label'], colnames=['predict'])
    print(cx)

    df1 = pd.DataFrame({'label':y_test_label,'predict':prediction})
    print(df1.shape)
    print(df1[:2])
    print(df1[(df1.label == 5) & (df1.predict == 3)])   # important

if __name__ == '__main__':
    main()
