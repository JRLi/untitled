from keras.models import Sequential
from keras.layers import Dense
from keras.utils.np_utils import to_categorical
from numpy import transpose
import numpy as np
import csv



with open("GSE68086.csv", "r", encoding = "utf-8") as f:
    rows = csv.reader(f)
    
    label_type_list = rows.__next__()[2:];
    
    types = list(set(label_type_list));
    type_index = {types[i]:i for i in range(len(types))}
    
    labels = [[type_index[lable_type]] for lable_type in label_type_list]
    labels = np.array(labels)
    labels = to_categorical(labels, len(type_index))
    
    rows.__next__();
    
    data = [row[2:] for row in rows]
    for row in data:
        for i in range(len(row)):
            row[i] = float(row[i])
    data = np.array(data)
    data = data / np.linalg.norm(data)
    data = transpose(data)



model = Sequential()
model.add(Dense(len(type_index), input_dim=len(data[0]), activation="softmax"))
 
model.compile(optimizer="rmsprop", loss="categorical_crossentropy", metrics=["accuracy"])

model.fit(data, labels, nb_epoch=10, batch_size=32)

score = model.predict_classes(data)

print()
print()
print("Type:", type_index)
print("All:", score)
print("Lung:", score[195:255])
print("HC:", score[121:181])