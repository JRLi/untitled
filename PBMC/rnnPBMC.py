from keras.models import Sequential, model_from_json
from keras.layers import Dense, Dropout
from PBMC.addonPBMCrnn import get_data, corr, save_model, load_model
import tensorflow as tf
import h5py
from sklearn.model_selection import train_test_split
tf.python.control_flow_ops = tf

output_dim1=64
output_dim2=32
batch_size=10
nb_epoch=1000
path = 'D:/Project/PBMC/logistic_in/'
file = 'GSE16129_19301_0_SAI_1_exacerbationAsthma'

data_train, labels_train, types_train, data_test, labels_test, types_test = get_data(path + file, 'Disease', 0.8)
print('type(data_train):', type(data_train))
print(data_train.shape)
print(types_train)
print(labels_train)
# print('data_train.iloc[1]:',data_train.iloc[2])
model = Sequential()

model.add(Dense(output_dim=output_dim1, input_dim=len(data_train[0]), activation="relu"))
model.add(Dropout(0.25))
model.add(Dense(output_dim2, activation="relu"))
model.add(Dropout(0.5))
model.add(Dense(len(types_train), activation="sigmoid"))

model.compile(loss="mse", optimizer="rmsprop")

model.fit(data_train, labels_train, nb_epoch=nb_epoch, batch_size=batch_size)

predict_types = model.predict_classes(data_test)
print("Types:", types_test)
print("True Types:", labels_test.argmax(1))
print("Predict Types:", predict_types)
print("Corr:", corr(labels_test.argmax(1).tolist(), predict_types.tolist()))


save_model(model, path)
model = load_model(path)

predict_types = model.predict_classes(data_test)
print("Types:", types_test)
print("True Types:", labels_test.argmax(1))
print("Predict Types:", predict_types)
print("Corr:", corr(labels_test.argmax(1).tolist(), predict_types.tolist()))

# select_sample = np.array([data[0]])
# score = model.predict_classes(select_sample)
