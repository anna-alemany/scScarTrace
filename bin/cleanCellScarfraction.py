import sys, os
import pandas as pd
from pandas.io.parsers import read_csv

try:
    df = read_csv(sys.argv[1], sep = '\t', index_col = 0)
    thr = float(sys.argv[2])
    outfile = sys.argv[3]
except:
    sys.exit('Please, give input df file to clean; threshold fraction to filter out; full name for output df file')

print df.shape

df = df[df>thr]
df = df.fillna(0)

scars = df.columns[:-1]

df = df.transpose()
df.loc[scars] = 100.*df.loc[scars]/df.loc[scars].sum()
df = df.transpose()
df['hclust'] = df['hclust'].astype(int)

df = df.iloc[:,(df.sum() > 0).values]

df.to_csv(outfile, sep = '\t')

print df.shape
