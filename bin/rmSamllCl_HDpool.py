import sys, os
from pandas.io.parsers import read_csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Colors import *

try:
    df = read_csv(sys.argv[1], sep = '\t', index_col = 0)
    out = sys.argv[2]
    pdfplot = sys.argv[3]
except:
    sys.exit('Please, give input HD_pool file (after running automat); output file; pdfplot option (y/n)')

print df.shape

wt = list(set(df[df['720M'] == 100]['hclust']))
if len(wt) > 1:
    sys.exit('more than one wt cluster!')
if len(wt) > 0:
    keepidx = list(df[df['hclust']==wt].index)
else:
    keepidx = []
    wt = [-1]
for cl in set(df['hclust']):
    if df[df['hclust'] == cl].shape[0] > 1 and cl != wt[0]:
        keepidx += list(df[df['hclust']==cl].index)

df = df.loc[keepidx]

print df.shape

df.loc[df[df['hclust']==wt[0]].index, 'hclust'] = -1

i = 0
dfnew = df.copy()
for cl in sorted(list(set(df['hclust']))):
        cells = df[df['hclust'] == cl].index
        dfnew.loc[cells, 'hclust'] = i
        i += 1
df = dfnew.copy()

print df.shape

scars = df.columns[:-1]
csel = []
for scar in scars:
    if df[scar].sum() > 0:
        csel += [scar]
csel += ['hclust']

df = df[csel]
print df.shape

df.to_csv(out, sep = '\t')
print i-1


if pdfplot=='y':
    fig = plt.figure(figsize=(15,5))
    bottom=np.zeros(len(df.index))
    for i,cigar in enumerate(df.columns[:-1]):
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

    fig.savefig(out + '_barplot.pdf',  additional_artist=art, bbox_inches='tight')

