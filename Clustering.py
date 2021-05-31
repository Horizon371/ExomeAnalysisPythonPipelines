import glob
import time

import numpy as np
import pandas as pd
from matplotlib import pyplot
from numpy import unique, where
from sklearn.cluster import KMeans

exAc1 = 'ExAC_ALL'
exAc2 = 'ExAC_AFR'

file_name = r"\clusters" + str(time.time()) + ".png"
save_path = r"D:\Projects\Licenta\Exome Analysis Project\saved clustering plots" + file_name

print(save_path)
path = r'D:\Projects\Licenta\Exome Analysis Project\saved exomes analysis'
all_files = glob.glob(path + "/*.csv")

X = []

for file_name in all_files:
    current_csv = pd.read_csv(file_name, index_col=None, na_values=['NA'], usecols=[exAc1, exAc2])

    filt1 = (current_csv[exAc1] != ".")
    filt2 = (current_csv[exAc2] != ".")

    ex_ac_all = current_csv.loc[filt1, [exAc1]]
    ex_ac_all[exAc1] = ex_ac_all[exAc1].astype(float)
    ex_ac_all_mean = ex_ac_all[exAc1].mean()

    ex_ac_eas = current_csv.loc[filt2, [exAc2]]
    ex_ac_eas[exAc2] = ex_ac_eas[exAc2].astype(float)
    ex_ac_eas_mean = ex_ac_eas[exAc2].mean()

    X.append([ex_ac_all_mean, ex_ac_eas_mean])

X = np.array(X)
Y = np.array([[1, 2], [1, 4], [1, 0], [10, 2], [10, 4], [10, 0]])

model = KMeans(n_clusters=2)
model.fit(X)
yhat = model.predict(X)
clusters = unique(yhat)
for cluster in clusters:
    row_ix = where(yhat == cluster)
    pyplot.scatter(X[row_ix, 0], X[row_ix, 1])
pyplot.savefig(save_path)
pyplot.plot()
