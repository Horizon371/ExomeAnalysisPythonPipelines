import glob
import sys

import pandas as pd

geneLocations = {}
geneName = sys.argv[1]

path = r'D:\Projects\Licenta\Exome Analysis Project\saved exomes analysis'
all_files = glob.glob(path + "/*.csv")

# now = round(time.time() * 1000)

print("Gene Frequency Pipeline started")

for file_name in all_files:
    truncatedDataFrame = pd.read_csv(file_name, index_col=None, na_values=['NA'], usecols=["GeneName", "CHROM", "POS"])

    filt = truncatedDataFrame['GeneName'] == geneName

    filteredDataFrame = truncatedDataFrame[filt]

    for index, row in filteredDataFrame.iterrows():
        if (row['POS'], row['CHROM']) not in geneLocations:
            geneLocations[(row['POS'], row['CHROM'])] = 1
        else:
            geneLocations[(row['POS'], row['CHROM'])] += 1


for key, value in geneLocations.items():
    print(key[0], key[1], value)


# then = round(time.time() * 1000)
# print(then-now)
