import pandas as pd
from pandas import DataFrame
from pandas import Series
from scipy import stats
import numpy as np
import os
from sklearn import svm

fc = 20.0
dataset = "Samples_combined_mRNA_TPM.txt"
#dataset = "Samples_combined_lncRNA_TPM.txt"
#dataset = "Samples_combined_circRNA_RPM.txt"
print(dataset)
df = pd.read_csv(os.path.join('D:/Project/lncRNA/exoRbase', dataset), sep="\t")

matchs = {"N": r"N\d{1,2}\Z", "CHD": r"CHD\d{1,2}\Z", "CRC": r"CRC\d{1,2}\Z", "HCC": r"HCC\d{1,2}\Z",
          "PAAD": r"PAAD\d{1,2}\Z", "WhB": r"WhB\d{1,2}\Z"}

ns = list(df.columns[df.columns.str.match(matchs["N"])])
chds = list(df.columns[df.columns.str.match(matchs["CHD"])])
crcs = list(df.columns[df.columns.str.match(matchs["CRC"])])
hccs = list(df.columns[df.columns.str.match(matchs["HCC"])])
paads = list(df.columns[df.columns.str.match(matchs["PAAD"])])
whbs = list(df.columns[df.columns.str.match(matchs["WhB"])])
print("n:{}, chd:{}, crc:{}, hcc:{}, paad:{}, whb:{}".format(len(ns), len(chds), len(crcs), len(hccs), len(paads),
                                                             len(whbs)))
label = list()
label_name = list()

ndf = DataFrame()

label_id = 0
label_name.append("ns")
for n in ns:
    ndf[len(label)] = df[n]
    label.append(label_id)
label_id = 1
label_name.append("chds")
for n in chds:
    ndf[len(label)] = df[n]
    label.append(label_id)
label_id = 2
label_name.append("crcs")
for n in crcs:
    ndf[len(label)] = df[n]
    label.append(label_id)
label_id = 3
label_name.append("hccs")
for n in hccs:
    ndf[len(label)] = df[n]
    label.append(label_id)
label_id = 4
label_name.append("paads")
for n in paads:
    ndf[len(label)] = df[n]
    label.append(label_id)
label_id = 5
label_name.append("whbs")
for n in whbs:
    ndf[len(label)] = df[n]
    label.append(label_id)

ndf.columns = label
clf = svm.SVC(kernel='linear')

type1 = 0
type2 = 3
msg = list()
runs = ndf[[type1, type2]].T
ids = list(set(runs.index))
print("{} vs {}".format(label_name[ids[0]], label_name[ids[1]]))
#print(runs)
#runs.to_csv(os.path.join('D:/Project/lncRNA/exoRbase', 'lncTest.csv'))
runs_miRNA = list()
tp, fp, fn, tn = 0, 0, 0, 0
print(runs.dtypes)
for i in range(len(runs)):
    print("Run {}".format(i), end="")
    run = runs[0:1]
    label = run.index[0]
    runs = runs[1:]     # change runs
    labels = runs.index
    train = runs.copy()
    test = run.copy()
    #print(label, labels)
    #print(np.array(labels))
    n = train.T[type1].replace(0.0, 0.01).mean(axis=1)  # FC
    t = train.T[type2].replace(0.0, 0.01).mean(axis=1)  # FC
    train = train.T[t > (fc * n)].T  # FC
    test = test.T[t > (fc * n)].T  # FC
    #     pvalues = train.apply(lambda x: stats.ttest_ind(x[0], x[1])[1]) # P Value
    #     train = train.T[pvalues < pvalue].T # P Value
    #     test = test.T[pvalues < pvalue].T # P Value
    #print(" #miRNA {}".format(len(train.columns)))
    runs_miRNA.append(len(train.columns))

    model = clf.fit(np.array(train), np.array(labels))
    predict = model.predict(np.array(test))[0]
    #print(label, predict)
    if label != 0 and predict != 0:
        tp += 1
    elif label == 0 and predict != 0:
        fp += 1
    elif label != 0 and predict == 0:
        fn += 1
    elif label == 0 and predict == 0:
        tn += 1

    runs = runs.append(run)     # recovery data frame

    p = 0
    r = 0
    try:
        p = (tp / (tp + fp))
    except:
        pass
    try:
        r = (tp / (tp + fn))
    except:
        pass
    #print(", TP {}, FP {}, TN {}, FN {}, precision {:.3f}, recall {:.3f}".format(tp, fp, tn, fn, p, r))
print()
precision = tp / (tp + fp)
recall = tp / (tp + fn)
msg.append(
    "{} vs {}: precision {:.3f} and recall {:.3f} FC {} #features: {:.3f}".format(label_name[ids[0]], label_name[ids[1]],
                                                                              precision, recall, fc,
                                                                              np.mean(runs_miRNA)))
print("\n".join(msg))
