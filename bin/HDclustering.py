import sys, os
from pandas.io.parsers import read_csv
import pandas as pd
import sklearn.cluster

try:
    df = read_csv(sys.argv[1], sep = '\t', index_col = 0)
    outfile = sys.argv[2]
except:
    sys.exit('Please, give path to _df.txt file (full name); root for output file')


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

sys.exit()
centroids = centroidsl5
hclust = sklearn.cluster.AgglomerativeClustering(n_clusters = 5)
hclust.fit(centroids[centroids.columns[:-1]])
hdf = pd.DataFrame({'cell': centroids.index, 'cl': hclust.labels_})
hdf = hdf.sort_values(by='cl')
centroids = centroids.loc[hdf['cell']]
centroids.index.name = 'cellID'
centroids.to_csv(outfile + '_HDcentroids-l5.txt', sep = '\t')
idxsel = [idx for idx in dfnew.index if dfnew.loc[idx, 'hclust'] in set(centroids['hclust'])]
dfnew.loc[idxsel].to_csv(outfile + '_HDl5.txt', sep = '\t')

centroids = centroidsg5
hclust = sklearn.cluster.AgglomerativeClustering(n_clusters = 5)
hclust.fit(centroids[centroids.columns[:-1]])
hdf = pd.DataFrame({'cell': centroids.index, 'cl': hclust.labels_})
hdf = hdf.sort_values(by='cl')
centroids = centroids.loc[hdf['cell']]
centroids.index.name = 'cellID'
centroids.to_csv(outfile + '_HDcentroids-g5.txt', sep = '\t')
idxsel = [idx for idx in dfnew.index if dfnew.loc[idx, 'hclust'] in set(centroids['hclust'])]
dfnew.loc[idxsel].to_csv(outfile + '_HDg5.txt', sep = '\t')

