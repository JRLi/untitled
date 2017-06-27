#!/usr/bin/env python3
from keras.datasets import cifar10
from keras.utils import np_utils
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation, Flatten
from keras.layers import Conv2D, MaxPooling2D, ZeroPadding2D
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
np.random.seed(10)


def plot_image_label_prediction(images, labels, label_dict, prediction, idx, num=10):
    fig=plt.gcf()
    fig.set_size_inches(12, 14)
    num = 25 if num > 25 else num   # set sub-graph number limitation
    for i in range(0, num):
        ax = plt.subplot(5, 5, 1 + i)   # build 5 x 5 sub-graph
        ax.imshow(images[idx], cmap='binary')   # plot sub-graph of idx
        title = 'fig{}, {}'.format(idx, label_dict[labels[idx][0]])     # set title of each sub-graph
        if len(prediction) > 0:     # prediction can be an empty list, indict no prediction
            title += '=>' + label_dict[prediction[idx]] # if import prediction result, set title with prediction result
        ax.set_title(title, fontsize=10)    # set title size
        ax.set_xticks([])   # remove the tick of x axis of sub-graph
        ax.set_yticks([])   # remove the tick of y axis of sub-graph
        idx += 1    # next index
    plt.show()


def show_predicted_probability(y, prediction, x_img, predicted_probability, label_dict, i):
    print('label:', label_dict[y[i][0]], 'predict:', label_dict[prediction[i]])
    plt.figure(figsize=(2, 2))
    plt.imshow(np.reshape(x_img[i], (32, 32, 3)))
    plt.show()
    for j in range(len(predicted_probability[i])):
        print(label_dict[j] + 'Probability:%1.9f'%(predicted_probability[i][j]))


def show_train_history(train_history, train, validation):
    plt.plot(train_history.history[train])
    plt.plot(train_history.history[validation])
    plt.title('Train History')
    plt.ylabel(train)
    plt.xlabel('Epoch')
    plt.legend(['train', 'validation'], loc='upper left')
    plt.show()


def cnn(X, y):
    model =Sequential()
    model.add(Conv2D(filters=32, kernel_size=(3, 3), input_shape=X[0].shape, activation='relu', padding='same'))
    model.add(Dropout(0.25))
    model.add(MaxPooling2D(pool_size=(2, 2)))
    model.add(Conv2D(filters=64, kernel_size=(3, 3), activation='relu', padding='same'))
    model.add(Dropout(0.25))
    model.add(MaxPooling2D(pool_size=(2, 2)))
    model.add(Flatten())
    model.add(Dropout(0.25))
    model.add(Dense(units=1024, activation='relu'))
    model.add(Dropout(0.25))
    model.add(Dense(units=y.shape[1], activation='softmax'))
    print(model.summary())
    model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
    train_history = model.fit(X, y, validation_split=0.2, epochs=10, batch_size=128, verbose=1)
    return train_history, model

def main():
    (X_img_train, y_label_train), (X_img_test, y_label_test) = cifar10.load_data()
    print(X_img_train.shape, X_img_test.shape, y_label_train.shape, y_label_test.shape)
    label_dict = {0: 'airplane', 1: 'automobile', 2: 'bird', 3: 'cat', 4: 'deer',
                  5: 'dog', 6: 'frog', 7: 'horse', 8: 'ship', 9: 'truck'}
    plot_image_label_prediction(X_img_train, y_label_train, label_dict, [], 15, 20) # shape of y isn't same with mnist
    # normalize X value
    X_img_train_normalize = X_img_train.astype('float32') / 255.0
    X_img_test_normalize = X_img_test.astype('float32') / 255.0
    y_label_train_OneHot = np_utils.to_categorical(y_label_train)
    y_label_test_OneHot = np_utils.to_categorical(y_label_test)
    print(X_img_train_normalize[0].shape, y_label_train_OneHot.shape)
    train_history, model = cnn(X_img_train_normalize, y_label_train_OneHot)
    show_train_history(train_history, 'acc', 'val_acc')
    show_train_history(train_history, 'loss', 'val_loss')

    scores = model.evaluate(X_img_test_normalize, y_label_test_OneHot, verbose=0)
    print(scores)
    print('accuracy:', scores[1])

    prediction = model.predict_classes(X_img_test_normalize)
    plot_image_label_prediction(X_img_test, y_label_test, label_dict, prediction, 0, 20)

    predicted_probability = model.predict(X_img_test_normalize)
    show_predicted_probability(y_label_test, prediction, X_img_test, predicted_probability, label_dict, 0)
    show_predicted_probability(y_label_test, prediction, X_img_test, predicted_probability, label_dict, 3)

    print(prediction.shape)     # 1d array
    print(y_label_test.shape)   # 2d array
    y_label_test_1d = y_label_test.reshape(-1)    # convert to 1d array
    print(label_dict)
    df_cm = pd.crosstab(y_label_test_1d, prediction, rownames=['label'], colnames=['predict'])
    print(df_cm)
    model.save_weights('E://StringTemp/cifarCnnModel.h5')

if __name__ == '__main__':
    main()
