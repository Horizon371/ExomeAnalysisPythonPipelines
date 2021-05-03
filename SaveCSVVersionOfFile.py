import sys

import pandas as pd

fileName = sys.argv[1]
xlsBaseDir = r'D:\Projects\Licenta\Exome Analysis Project\saved exomes\\'
csvBaseDir = r'D:\Projects\Licenta\Exome Analysis Project\saved exomes analysis\\'
df = pd.read_excel(xlsBaseDir + fileName + '.xls')
df.to_csv(csvBaseDir + fileName + '.csv')