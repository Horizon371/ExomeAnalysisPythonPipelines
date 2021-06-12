import sys

import pandas as pd

genes = ['AIP', 'ALK', 'APC', 'ATM', 'AXIN2', 'BAP1', 'BARD1', 'BLM', 'BMPR1A', 'BRCA1', 'BRCA2', 'BRIP1', 'CASR',
         'CDC73', 'CDH1', 'CDK4', 'CDKN1B', 'CDKN1C', 'CDKN2A', 'CEBPA', 'CHEK2', 'CTNNA1', 'DICER1', 'DIS3L2', 'EGFR',
         'EPCAM', 'FH', 'FLCN', 'GATA2', 'GPC3', 'GREM1', 'HOXB13', 'HRAS', 'KIT', 'MAX', 'MEN1', 'MET', 'MITF', 'MLH1',
         'MSH2', 'MSH3', 'MSH6', 'MUTYH', 'NBN', 'NF1', 'NF2', 'NTHL1', 'PALB2', 'PDGFRA', 'PHOX2B', 'PMS2', 'POLD1',
         'POLE', 'POT1', 'PRKAR1A', 'PTCH1', 'PTEN', 'RAD50', 'RAD51C', 'RAD51D', 'RB1', 'RECQL4', 'RET', 'RUNX1',
         'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCE1', 'STK11', 'SUFU', 'TERC',
         'TERT', 'TMEM127', 'TP53', 'TSC1', 'TSC2', 'VHL', 'WRN', 'WT1']

fileName = sys.argv[1]
baseDir = r'D:\Projects\Licenta\Exome Analysis Project\saved exomes analysis\\'
df = pd.read_csv(baseDir + fileName + ".csv")

filt = (df['GeneName'].isin(genes)) & (df['OMIM'] != ".") & (df['Func'] != "intergenic") & (
        df['Func'] != "intronic") & (
               df['Func'] != "UTR5") & (
               df['Func'] != "UTR3") & (
               df['Func'] != "downstream") & (
               df['Func'] != "upstream")

with pd.option_context('display.max_rows', None, 'display.max_columns', None):
    print(df.loc[filt, ['GeneName', 'CHROM', 'POS']])
