import sys

import pandas as pd

fileName = sys.argv[1]
baseDir = r'D:\Projects\Licenta\Exome Analysis Project\saved exomes analysis\\'
df = pd.read_csv(baseDir + fileName + ".csv")
filt = (df['Ex52'].str.slice(stop=3) == '1/1')
percentageTable = filt.value_counts(normalize=True)
print(percentageTable[1])

