from keras.models import Sequential, model_from_json
from keras.layers import Dense, Dropout
from PBMC.addonPBMCrnn import get_data, corr, save_model, load_model, score, generate_results, roc_plot
import tensorflow as tf
import os
tf.python.control_flow_ops = tf
path = 'D:/Project/PBMC/logistic_in/'
out_path = 'D:/Project/PBMC/logistic_out/'
file = 'GSE16129_19301_0_SAI_1_exacerbationAsthma'
fpr = dict()
tpr = dict()
roc_auc = dict()

output_dim1=64
output_dim2=32
batch_size=10
nb_epoch=1000

# if the data format is not standard, must assign the mode parameter to not 1.
data_train, labels_train, types_train, data_test, labels_test, types_test = get_data(path + file, 'Disease', 0.8)
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

roc_plot(fpr, tpr, roc_auc)
