from keras.models import model_from_json
from keras.utils.np_utils import to_categorical
from numpy import transpose
import numpy as np
import scipy.stats as stats
import csv


def get_data(path):
    with open(path, "r", encoding="utf-8") as f:
        rows = csv.reader(f)

        label_type_list = rows.__next__()[2:]
        print('rows.__next__()[2:]: ', label_type_list)  # debug
        types = list(set(label_type_list))
        types.sort()
        print('types after sort:', types)    # debug
        types = {types[i]: i for i in range(len(types))}
        print('types_dict:', types)
        labels = [[types[lable_type]] for lable_type in label_type_list]
        labels = np.array(labels)
        labels = to_categorical(labels, len(types))

        rows.__next__()

        data = [row[2:] for row in rows]
        for row in data:
            for i in range(len(row)):
                row[i] = float(row[i])
        data = np.array(data)
        data = data / np.linalg.norm(data)
        data = transpose(data)

    return data, labels, types


def corr(a, b):
    return stats.pearsonr(a, b)[0]


def save_model(model):
    json_string = model.to_json()
    with open("model.json", "w") as output:
        output.writelines(json_string)
    model.save_weights("weights.h5")


def load_model():
    with open("model.json", "r") as input:
        model = model_from_json(input.readline())
    model.load_weights("weights.h5")
    return model
