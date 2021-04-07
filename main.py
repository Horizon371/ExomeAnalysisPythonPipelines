import numpy as np
import pandas as pd

data_file = pd.read_excel(r'D:\Projects\Licenta\Exome Analysis Project\saved exomes\EXOM1.xls')
print(data_file['GeneName'])