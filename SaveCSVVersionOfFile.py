import sys

import pandas as pd

fileName = sys.argv[1]
xlsBaseDir = r'D:\Projects\Licenta\Exome Analysis Project\saved exomes\\'
csvBaseDir = r'D:\Projects\Licenta\Exome Analysis Project\saved exomes analysis\\'
print(xlsBaseDir + fileName + '.xls')
df = pd.read_excel(xlsBaseDir + fileName + '.xls')
df.columns = df.columns.str.replace('Ex[0-9]+', 'Ex52')
df.columns = df.columns.str.replace('EX[0-9]+', 'Ex52')
df.columns = df.columns.str.replace('ex[0-9]+', 'Ex52')
df.to_csv(csvBaseDir + fileName + '.csv')