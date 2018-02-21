import sys, os
from pandas.io.parsers import read_csv
import numpy as np
import pandas as pd
import sklearn.cluster
import matplotlib.pyplot as plt
from Colors import *

try:
    df = read_csv(sys.argv[1], sep = '\t', index_col = 0)
    outfile = sys.argv[2]
    pdfplot = sys.argv[3]
except:
    sys.exit('Please, give path to _df.txt file (full name); root for output file; produce pdf plot (y/n)')

hcl = df['hclust']
del df['hclust']

clones = {}
for cell in df.index:
#    scars = '-'.join(df.columns[df.loc[cell] > 0])
    scars = '-'.join(df.columns[df.loc[cell] > 3.5])
    if scars in clones:
        clones[scars].append(cell)
    else:
        clones[scars] = [cell]

dfnew = pd.DataFrame()
i = 0
for clone in clones:
    pdf = df.loc[clones[clone]]
    pdf['hclust'] = i
    dfnew = pd.concat([dfnew, pdf])
    i += 1

dfnew.to_csv(outfile + '_HD.txt', sep = '\t')

f = open(outfile + "_HD.gpl", 'w')
print >> f, 'set style data histogram'
print >> f, 'set style histogram rowstacked'
print >> f, 'set style fill solid'
print >> f, 'set key below'
print >> f, 'unset xtic'
print >> f, 'set ytics 0,20,100'
print >> f, 'pl for [i=2:'+str(dfnew.shape[1]-1)+'] "'+outfile+'_HD.txt" us i:xtic(1) ti col'
f.close()


centroids = pd.DataFrame(columns = dfnew.columns)
centroidsl5 = pd.DataFrame(columns = dfnew.columns)
centroidsg5 = pd.DataFrame(columns = dfnew.columns)
for clone in set(dfnew['hclust']):
    centroids.loc[clone] = dfnew[dfnew['hclust']==clone].mean()
    if dfnew[dfnew['hclust'] == clone].shape[0] < 5:
        centroidsl5.loc[clone] = dfnew[dfnew['hclust'] == clone].mean()
    else:
        centroidsg5.loc[clone] = dfnew[dfnew['hclust'] == clone].mean()

hclust = sklearn.cluster.AgglomerativeClustering(n_clusters = 5)
hclust.fit(centroids[centroids.columns[:-1]])
hdf = pd.DataFrame({'cell': centroids.index, 'cl': hclust.labels_})
hdf = hdf.sort_values(by='cl')
centroids = centroids.loc[hdf['cell']]
centroids.index.name = 'cellID'
centroids.to_csv(outfile + '_HDcentroids.txt', sep = '\t')


if pdfplot=='y':
    fig = plt.figure(figsize=(15,5))
    bottom=np.zeros(len(dfnew.index))
    for i,cigar in enumerate(dfnew.columns[:-1]):
        j = np.mod(i,len(colors))
        plt.bar(range(len(dfnew.index)), dfnew[cigar], bottom = bottom, width = 1, color=colors[j])
        bottom += dfnew[cigar]

    plt.ylim(0,100)
    plt.xlim(0,len(dfnew.index))
    plt.ylabel('scar %')
    plt.xlabel('cells')
    
    art = []
    lgd = plt.legend(df.columns[:-1], loc = 9, bbox_to_anchor = (0.5, -0.1), ncol = 5)
    art.append(lgd)

    fig.savefig(outfile + '_HDbarplot.pdf',  additional_artist=art, bbox_inches='tight')
    
