import sys, os
import sklearn.cluster
from pandas.io.parsers import read_csv
import numpy as np
import pandas as pd
import scipy.spatial.distance as dist

try:
    indata = sys.argv[1]
    ncluster = int(sys.argv[2])
    outfile = sys.argv[3]
    method = sys.argv[4]
except:
    sys.exit('Please, give input table, number of clusters; root for output file; method (hcl/acl)')

df = read_csv(indata, sep = '\t', index_col = 0)
df.columns.name = 'cellid'
if '720M' in df.index:
    df = df.loc[['720M'] + [idx for idx in df.index if idx != '720M']]

if method == 'hcl':  #### Before cleaning ####
    hclust = sklearn.cluster.KMeans(n_clusters = ncluster, random_state = 20)
    hclust.fit(df.transpose())
elif method == 'acl': #### After cleaning ####
    hclust = sklearn.cluster.AgglomerativeClustering(n_clusters = ncluster)
    hclust.fit(df.transpose())

hclustdf = pd.DataFrame({'hclust': hclust.labels_}, index = df.columns)
hclustdf = hclustdf.sort_values(by = 'hclust')
hclustdf.to_csv(outfile + '_clust.txt', sep = '\t', header = None)

centroiddf = pd.DataFrame(index = df.index)
centroiddf.columns.name = 'centroids'
distCM = pd.DataFrame(columns = ['dmeanCM', 'dvarCM', 'numCell'])
for i in range(ncluster):
    cells = hclustdf[hclustdf['hclust'] == i].index
    rdf = df[cells]
    centroid = rdf.mean(axis=1).sort_values(ascending=False)
    centroiddf[i] = centroid
    d = []
    for cell in cells:
        d.append(dist.euclidean(df[cell], centroid))
    distCM.loc[i] = [np.array(d).mean(), np.array(d).var(), len(cells)]
centroiddf.transpose().to_csv(outfile + '_centroid.txt', sep = '\t')
distCM.to_csv(outfile + '_distCM.txt', sep = '\t')

intSize0 = centroiddf.values.tolist()
intSize = []
for l in intSize0:
    intSize += l
binnum = 90
h = np.histogram(intSize, bins = binnum, density = True)
histo = pd.DataFrame({'x': [ 0.5*sum(h[1][i:i+2]) for i in range(binnum)], 'y': h[0]})
histo.to_csv(outfile + '_scarSize.txt', sep = '\t', index = None)

df[hclustdf.index].corr().to_csv(outfile + '_corr.txt', sep = '\t')

df = df[hclustdf.index].transpose()
df['hclust'] = hclustdf
df.to_csv(outfile + '_df.txt', sep = '\t')

f = open(outfile + '.gpl', 'w')
print >> f, 'set style data histogram'
print >> f, 'set style histogram rowstacked'
print >> f, 'set style fill solid'
print >> f, 'unset xtic'
print >> f, 'set key out'
print >> f, 'set ytics 0,20,100'
print >> f, 'l "var/gnuplot-extendcolor.gpl"'
print >> f, 'pl for [i=2:' + str(df.shape[1]) + '] "' + outfile + '_df.txt" us i:xtic(1) ti col'
print >> f, 'rep for [i=0:' + str(ncluster-1) + '] "' + outfile + '_clust.txt" us ($2==i?10:1/0) noti'
f.close()



