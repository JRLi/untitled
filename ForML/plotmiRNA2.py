import seaborn as sns
import matplotlib.pyplot as plt
import csv
import re
import random as r
import numpy as np
from scipy import stats
from scipy.integrate import quad
import warnings
warnings.filterwarnings("ignore")
match_map = {"Prostate": [r"TCGA-PRAD", r"PC\d{1,2}\Z"], "Colon": [r"TCGA-COAD", r"\d{1}S\d{1,2}"], "Healthy": [r"N\d{1,2}\Z"]}
target_cancer = "Prostate"
target_miRNA = "miR-375"
target_theta = 0.2
#target_cancer = input("Cancer(enter Prostate or Colon): ")
#target_miRNA = input("miRNA(ex. miR-375): ")
#target_theta = eval(input("Theta(ex. 0.2): "))
with open("D:/Project/circ_miRNA/TCGA.csv", "r") as inp:
    csv_inp = csv.reader(inp)
    cancer, label, sample, *miRNAs = csv_inp
    for miRNA in miRNAs:
        if target_miRNA.lower() in miRNA[0].lower():
            tcga_exp = []
            tcga_miRNA = miRNA[0]
            matchs = match_map[target_cancer]
            for i in range(len(cancer)):
                for match in matchs:
                    pattern = re.compile(match)
                    match = pattern.match(cancer[i])
                    if match != None:
                        if label[i] == "Cancer":
                            tcga_exp.append(float(miRNA[i]))

geo_cancer_exps = list()
geo_normal_exps = list()
geo_miRNAs = list()
with open("D:/Project/circ_miRNA/GSE71008_Data_matrix.csv", "r") as inp:
    csv_inp = csv.reader(inp)
    title, *miRNAs = csv_inp
    for miRNA in miRNAs:
        if target_miRNA.lower() in miRNA[0].lower():
            exps = list()

            matchs = match_map[target_cancer]
            for i in range(len(title)):
                for match in matchs:
                    pattern = re.compile(match)
                    match = pattern.match(title[i])
                    if match != None:
                        exps.append(float(miRNA[i]))

            geo_cancer_exps.append(exps)

            exps = list()
            matchs = match_map["Healthy"]
            for i in range(len(title)):
                for match in matchs:
                    pattern = re.compile(match)
                    match = pattern.match(title[i])
                    if match != None:
                        exps.append(float(miRNA[i]))

            geo_normal_exps.append(exps)
            geo_miRNAs.append(miRNA[0])
print(geo_cancer_exps)

for i in range(len(geo_miRNAs)):
    print(geo_miRNAs[i])
    f, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(8, 9))
    ax1.set_title("{}\n{} cancer plasma (mean={:.2f}, sd={:.2f}, N={})".format(target_cancer, geo_miRNAs[i],
                                                                               np.mean(geo_cancer_exps[i]),
                                                                               np.std(geo_cancer_exps[i]),
                                                                               len(geo_cancer_exps[i])))
    #sns.set(style="white", palette="muted", color_codes=True)
    ax2.set_title("{} normal plasma (mean={:.2f}, sd={:.2f}, N={})".format(geo_miRNAs[i], np.mean(geo_normal_exps[i]),
                                                                           np.std(geo_normal_exps[i]),
                                                                           len(geo_normal_exps[i])))
    #sns.set(style="white", palette="muted", color_codes=True)
    ax3.set_title("{} tumor (mean={:.2f}, sd={:.2f}, N={})".format(tcga_miRNA, np.mean(tcga_exp), np.std(tcga_exp),
                                                                   len(tcga_exp)))
    sns.set(style="white", palette="muted", color_codes=True)
    sns.distplot(geo_cancer_exps[i], hist=False, color="b", ax=ax1)
    sns.distplot(geo_normal_exps[i], hist=False, color="g", ax=ax2)
    sns.distplot(tcga_exp, hist=False, color="r", ax=ax3)
    simu_exps = list()
    for exp in geo_cancer_exps[i]:
        rindex = r.randint(0, len(tcga_exp) - 1)
        simu_exps.append(((1 - target_theta) * exp + target_theta * tcga_exp[rindex]))
    ax4.set_title(
        "theta {} distribution (mean={:.2f}, sd={:.2f}, N={})".format(target_theta, np.mean(simu_exps), np.std(simu_exps),
                                                                      len(simu_exps)))
    #sns.set(style="white", palette="muted", color_codes=True)
    sns.distplot(simu_exps, hist=False, color="m", ax=ax4)

    overlap_list = list()
    for test_theta in range(1, 100):
        test_theta *= 0.01
        test_exps = list()
        for exp in geo_cancer_exps[i]:
            rindex = r.randint(0, len(tcga_exp) - 1)
            test_exps.append(((1 - test_theta) * exp + test_theta * tcga_exp[rindex]))

        a = np.array(simu_exps)
        b = np.array(test_exps)

        ker_a = stats.gaussian_kde(a)
        ker_b = stats.gaussian_kde(b)

        overlap = quad(lambda pt: min(ker_a(pt), ker_b(pt)), -np.inf, np.inf)

        overlap_list.append((overlap[0], test_theta, test_exps))
        print("theta {:.2f} overlap: {}".format(test_theta, overlap[0]))

    plt.tight_layout()

    f, axs = plt.subplots(3, 1, figsize=(8, 9))
    overlap_list.sort(reverse=True)
    for j in range(len(axs)):
        overlap, test_theta, test_exps = overlap_list[j]
        if j == 0:
            axs[j].set_title(
                "{} {} Top 3 Overlap\ntheta {} distribution {:.2f} overlap".format(target_cancer, geo_miRNAs[i],
                                                                                   test_theta, overlap))
        else:
            axs[j].set_title("theta {} distribution {:.2f} overlap".format(test_theta, overlap))
        axs[j].legend(ncol=2, loc="upper right", frameon=True)
        sns.set(style="white", palette="muted", color_codes=True)
        sns.distplot(simu_exps, hist=False, color="m", ax=axs[j], label="Original")
        sns.distplot(test_exps, hist=False, color="c", ax=axs[j], label="Test")

    plt.tight_layout()
plt.show()
