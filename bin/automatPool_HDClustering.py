import sys, os
import pandas as pd
import numpy as np
from pandas.io.parsers import read_csv
import scipy.spatial.distance as dist
import matplotlib.pyplot as plt
from Colors import *

try:
    df = read_csv(sys.argv[1], sep = '\t', index_col=0)
    outfile = sys.argv[2]
    pdfplot=sys.argv[3]
except:
    sys.exit('Please, give input HD.txt file; and name of output file; pdf plot option (y/n)')

alpha = read_csv('/Users/anna/Desktop/zebrafish/dynamics/mRNA/alphaValues.txt', sep = '\t', index_col = 0)

def findNeigh(clone, clonedf, obsdf):
    scclone = clonedf.columns[clonedf.loc[clone]>0]
    neigh = []
    for cl in clonedf.index:
        if cl != clone and obsdf.loc[cl,'n']>=2:
            sccl = clonedf.columns[clonedf.loc[cl]>0]
            a = True
            for sc in scclone:
                a *= sc in sccl
            if a:
                neigh.append(cl)
    return neigh
k0 = True
while k0:
    scars = df.columns[:-1]
    scarsnowt = [s for s in scars if s != '720M']

    clonedf = pd.DataFrame(df[df['hclust']==i][scars].mean() for i in set(df['hclust']))
    obsdf = pd.DataFrame(df[df['hclust']==i].shape[0] for i in set(df['hclust']))
    obsdf.columns = ['n']

    probscar = {'720M': 1}
    for scar in scars:
        if scar != '720M':
            if scar in alpha.index:
                p = alpha.loc[scar, 'a_mean'] if  alpha.loc[scar, 'a_mean'] > 1e-5 else 1e-5
            else:
                p = 1e-5
            probscar[scar] = p
    probscar = pd.DataFrame(probscar.values(), index = probscar.keys(), columns = ['prob'])

    probclone = pd.DataFrame(columns = ['prob'])
    for clone in clonedf.index:
        probclone.loc[clone] = (probscar.loc[clonedf.columns[clonedf.loc[clone] > 0]].prod()) * (0.06**(clonedf.loc[clone,scarsnowt]>0).sum())

    pmis = 0.9
    merged = 0
    multiple = []
    for c1 in clonedf.index:
        c2 = findNeigh(c1, clonedf, obsdf)
        if len(c2) == 1 and clonedf.loc[c1, '720M'] < 100.:
            c2 = c2[0]
            scardf = ((clonedf.loc[c2]>0).astype(int) - (clonedf.loc[c1]>0).astype(int))
#           if probclone.loc[c1,'prob'] < pmis**(scardf*obsdf.loc[c1,'n']):
            if (scardf == -1).sum() == 0:
                merged += 1
                df.loc[df[df['hclust']==c1].index,'hclust'] = c2
                print c1, ' merged with ', c2
        elif len(c2) > 1 and clonedf.loc[c1, '720M'] < 100:
            print c1, 'has many choices:', c2

    print '--'
    df = df.sort_values(by='hclust')
    n = 0
    for c in set(df['hclust']):
        df.loc[df[df['hclust']==c].index, 'hclust'] = n
        n += 1


    if merged == 0:
        k0 = False


df.to_csv(outfile + '.txt', sep = '\t')

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

    fig.savefig(outfile + '_barplot.pdf',  additional_artist=art, bbox_inches='tight')

