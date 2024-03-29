import glob

import numpy as np
import pandas as pd
from matplotlib import pyplot
from numpy import unique, where

from KMeansImpl import KMeansImpl

genes = ['AIP', 'ALK', 'APC', 'ATM', 'AXIN2', 'BAP1', 'BARD1', 'BLM', 'BMPR1A', 'BRCA1', 'BRCA2', 'BRIP1', 'CASR',
         'CDC73', 'CDH1', 'CDK4', 'CDKN1B', 'CDKN1C', 'CDKN2A', 'CEBPA', 'CHEK2', 'CTNNA1', 'DICER1', 'DIS3L2', 'EGFR',
         'EPCAM', 'FH', 'FLCN', 'GATA2', 'GPC3', 'GREM1', 'HOXB13', 'HRAS', 'KIT', 'MAX', 'MEN1', 'MET', 'MITF', 'MLH1',
         'MSH2', 'MSH3', 'MSH6', 'MUTYH', 'NBN', 'NF1', 'NF2', 'NTHL1', 'PALB2', 'PDGFRA', 'PHOX2B', 'PMS2', 'POLD1',
         'POLE', 'POT1', 'PRKAR1A', 'PTCH1', 'PTEN', 'RAD50', 'RAD51C', 'RAD51D', 'RB1', 'RECQL4', 'RET', 'RUNX1',
         'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCE1', 'STK11', 'SUFU', 'TERC',
         'TERT', 'TMEM127', 'TP53', 'TSC1', 'TSC2', 'VHL', 'WRN', 'WT1']

exAc = 'ExAC_ALL'


file_name = r"\clusters" + str(1) + ".png"
save_path = r"D:\Projects\Licenta\Exome Analysis Project\saved clustering plots" + file_name

print(save_path)
path = r'D:\Projects\Licenta\Exome Analysis Project\saved exomes analysis'
all_files = glob.glob(path + "/*.csv")


def create_data(X):
    for file in all_files:
        try:
            current_csv = pd.read_csv(file, index_col=None, na_values=['NA'], usecols=[exAc, 'GeneName'])
            number_of_total_genes = current_csv.shape[0]
            filt = (current_csv['GeneName'].isin(genes))
            filtered_data_frame = current_csv.loc[filt, [exAc]]
            print(get_name_of_file(file))
            ex_ac_all_mean = get_mean_for(filtered_data_frame, exAc)
            number_of_cancer_genes = filtered_data_frame.shape[0]
            X.append([ex_ac_all_mean, number_of_cancer_genes/number_of_total_genes*100])
        except BaseException:
            pass
    return X


def get_file_names(file_names):
    for file in all_files:
        file_names.append(get_name_of_file(file))
    return file_names


def get_name_of_file(file):
    start = file.rfind('\\')
    name = file[start+1:]
    end = name.find(".")
    return name[:end]


def get_mean_for(current_csv, exAc):
    filt = (current_csv[exAc] != ".")
    ex_ac = current_csv.loc[filt, [exAc]]
    ex_ac[exAc] = ex_ac[exAc].astype(float)
    # print(ex_ac.shape[0], ex_ac[exAc].median())
    return ex_ac[exAc].median()


def write_to_log_file(cluster):
    f = open('clustering_logs.txt', 'a')
    f.write(cluster + "\n")
    f.close()


X = np.array(create_data([]))
file_names = np.array(get_file_names([]))


model = KMeansImpl(k=2)
model.fit(X)
yhat = model.predict(X)
clusters = unique(yhat)
for cluster in clusters:
    row_ix = where(yhat == cluster)
    write_to_log_file(str(cluster) + " : " + str(file_names[row_ix]))
    pyplot.scatter(X[row_ix, 0], X[row_ix, 1])

# pyplot.scatter(X[:, 0], X[:, 1])

pyplot.xlabel("ExAc_All")
pyplot.ylabel("Pathogenic variant percentage")
pyplot.plot()
pyplot.savefig(save_path)
print("Clustering done")
# pyplot.show()
