import sys, os
import pandas as pd
from pandas.io.parsers import read_csv

try:
    inputfile=sys.argv[1]
    cellthr=int(sys.argv[2])
    scarthr=int(sys.argv[3])
    out=sys.argv[4]
except:
    sys.exit('Please, give input, cell threshold, scar threshold, output root')

df = read_csv(inputfile, sep = '\t', index_col=0)
print df.sum().sum(), ' total reads'
print df.shape, 'initial shape'
f1 = df[df.columns[df.sum() > cellthr]]

fdf = f1.loc[f1.index[f1.sum(axis=1) > scarthr]]
fdf = fdf[fdf.columns[fdf.sum()>0]]
fdf = fdf.loc[fdf.sum(axis=1).sort_values(ascending=False).index]

ndf = 100.*fdf/fdf.sum()

#ndf = ndf.loc[ndf.sum(axis=1).sort_values(ascending=False).index]
fdf.to_csv(out+'-filter.txt', sep = '\t')
ndf.to_csv(out+'-norm-filter.txt', sep = '\t')

print fdf.shape, 'shape after filtering'
