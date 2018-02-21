import sys, os
import pandas as pd
import numpy as np
from pandas.io.parsers import read_csv
from Colors import *
import matplotlib.pyplot as plt

try:
    df = read_csv(sys.argv[1], sep = '\t', index_col = 0)
    thr = float(sys.argv[2])
    outfile = sys.argv[3]
    pdfplot=sys.argv[4]
except:
    sys.exit('Please, give input df file to clean; threshold fraction to filter out; root name for output df file; generation of pdfplot (y/n)')

print df.shape

df = df[df>thr]
df = df.fillna(0)

scars = df.columns[:-1]

df = df.transpose()
df.loc[scars] = 100.*df.loc[scars]/df.loc[scars].sum()
df = df.transpose()
df['hclust'] = df['hclust'].astype(int)

df = df.iloc[:,(df.sum() > 0).values]

df.to_csv(outfile + '_df.txt', sep = '\t')

print df.shape

if pdfplot=='y':
    fig = plt.figure(figsize=(15,5))
    bottom=np.zeros(len(df.index))
    for i, cigar in enumerate(df.columns[:-1]):
        j = np.mod(i,len(colors))
        plt.bar(range(len(df.index)), df[cigar], bottom = bottom, width = 1, color=colors[j])
        bottom += df[cigar]

    plt.ylim(0,100)
    plt.xlim(0,len(df.index))
    plt.ylabel('scar %')
    plt.xlabel('cells')

    art = []
    lgd = plt.legend(df.columns[:-1], loc = 9, bbox_to_anchor = (0.5, -0.1), ncol = 5)
    art.append(lgd)

    fig.savefig(outfile + '_barplot.pdf',  additional_artist=art, bbox_inches='tight')

