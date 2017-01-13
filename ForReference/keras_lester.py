# %%
# from keras.datasets import mnist
from keras.models import Sequential
from keras.layers.core import Dense, Dropout, Activation, Flatten
from keras.layers.convolutional import Convolution2D, MaxPooling2D
from keras.utils import np_utils

# to plot import matplotlib
import matplotlib
import matplotlib.pyplot as plt

# others
import h5py as h5
import numpy as np
import collections as col

# %%
batch_size = 128
nb_epoch = 2

img_rows, img_cols = 38, 4
nb_filter = 8
filter_len = 4
nb_pool = 1

# %%
# define data
f_train = h5.File('D:/Project/PBMC/tools/machine_learning-master/ScrWT_train.h5', 'r')
f_test = h5.File('D:/Project/PBMC/tools/machine_learning-master/ScrWT_test.h5', 'r')
data_train = col.OrderedDict()

data_test = col.OrderedDict()
for k, v in f_train['data'].items():
    data_train[k] = v
for k, v in f_test['data'].items():
    data_test[k] = v

X_train = np.array(data_train['s_x'])
X_train = X_train.reshape(X_train.shape[0] ,1 ,X_train.shape[1] ,4)
Y_train = np.array(data_train['c0_y'])
# indices = Y_train < 0.7
# Y_train[indices] = 0
# indices = Y_train > 0
# Y_train[indices] = 1

X_test = np.array(data_test['s_x'])
X_test = X_test.reshape(X_test.shape[0] ,1 ,X_test.shape[1] ,4)
Y_test = np.array(data_test['c0_y'])
# indices = Y_test < 0.7
# Y_test[indices] = 0
# indices = Y_test > 0
# Y_test[indices] = 1


print('seq:', data_train['s_x'][0].shape)
print('target:', data_train['c0_y'][0])
print('X_train::', X_train[0])
print('Y_train:', Y_train[1 ,:])

# %%
# define the model
model = Sequential()

# add first convolution layer
model.add(Convolution2D(nb_filter, filter_len, filter_len,
                        border_mode='same',
                        input_shape=(1, img_rows, img_cols)))
convout = Activation('relu')
model.add(convout)

# add pooling layer
model.add(MaxPooling2D(pool_size=(nb_pool, nb_pool)))
model.add(Dropout(0.25)) # regulizer

model.add(Flatten())
model.add(Dense(128)) # add 128 neurons
model.add(Activation('relu'))
model.add(Dropout(0.5))
model.add(Dense(1))
model.add(Activation('sigmoid'))

model.compile(loss='mse', optimizer='rmsprop')
# model.compile(loss='binary_crossentropy', optimizer='adam')


# %%
# train the model
model.fit(X_train, Y_train, batch_size=batch_size, nb_epoch=nb_epoch,
          show_accuracy=False, verbose=1, validation_data=(X_test, Y_test))

# model.fit(X_train, Y_train, batch_size=batch_size, nb_epoch=nb_epoch,
#          show_accuracy=True, verbose=1, validation_split=0.2)

# %%
# evaluation
import sklearn.metrics as skm

def auc(y, z):
    if len(y) == 0 or len(np.unique(y)) < 2:
        return np.nan
    return skm.roc_auc_score(y, z)

def mse(y, z):
    return np.mean((y - z ) **2)

def rmse(y, z):
    return np.sqrt(mse(y, z))

def rrmse(y, z):
    return 1 - rmse(y, z)


z = model.predict(X_test, batch_size)
rrmse(Y_test, z)

