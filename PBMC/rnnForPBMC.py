from keras.models import Sequential, model_from_json
from keras.layers import Dense, Dropout
from PBMC.addonForRNN import get_data, corr, save_model, load_model
import tensorflow as tf
import h5py
from sklearn.model_selection import train_test_split
tf.python.control_flow_ops = tf

output_dim1=1024
output_dim2=512
batch_size=200
nb_epoch=1000
path = 'D:/Project/python_learning/'
training = 'Training.csv'
testing = 'Test.csv'

data, labels, types = get_data(path + training)  # 1: feature list in column 1, 2 or other: feature list in row 1
print(data.shape)
model = Sequential()
print(types)
model.add(Dense(output_dim=output_dim1, input_dim=len(data[0]), activation="relu"))
model.add(Dropout(0.25))
model.add(Dense(output_dim2, activation="relu"))
model.add(Dropout(0.5))
model.add(Dense(len(types), activation="sigmoid"))

model.compile(loss="mse", optimizer="rmsprop")

model.fit(data, labels, nb_epoch=nb_epoch, batch_size=batch_size)


data2, labels2, types2 = get_data(path + testing)
predict_types = model.predict_classes(data2)
print("Types:", types2)
print("True Types:", labels2.argmax(1))
print("Predict Types:", predict_types)
print("Corr:", corr(labels2.argmax(1).tolist(), predict_types.tolist()))


save_model(model, path)
model = load_model(path)

predict_types = model.predict_classes(data2)
print("Types:", types2)
print("True Types:", labels2.argmax(1))
print("Predict Types:", predict_types)
print("Corr:", corr(labels2.argmax(1).tolist(), predict_types.tolist()))

# select_sample = np.array([data[0]])
# score = model.predict_classes(select_sample)
