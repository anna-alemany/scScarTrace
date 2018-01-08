import sys, os
import itertools as it
from pandas.io.parsers import read_csv
import numpy as np

try:
    inputdffile = sys.argv[1]
    path2cigars = sys.argv[2]
    outputfile = sys.argv[3]
except:
    sys.exit('Please, give (1) input dataframe, filtered and normalized; (2) path to cigar seq-content; (3) root name for output file')

df = read_csv(inputdffile, sep = '\t', index_col = 0)

scarthreshold = 0.5
hdthreshold = 5

def HammingDistance(s1, s2):
    if s1.shape[0] != s2.shape[0]:
        xmin = min([s1.shape[0], s2.shape[0]])
        s1 = s1.loc[range(xmin)]
        s2 = s2.loc[range(xmin)]
#        sys.exit('Problemo in cigar-sequence shapes')
    hd1 = 0
    hd2 = 0
    for i in range(s1.shape[0]):
        if s1.loc[i, 'seq'] != s2.loc[i, 'seq']:
            hd1 += 1
        hd2 += np.sqrt(sum([(s1.loc[i, nt]-s2.loc[i,nt])**2 for nt in ['A', 'C', 'T', 'G']]))
    return hd1, hd2

minr = df.sum(axis=1).min() - 0.1

f = open(outputfile + '.log', 'w')

scarsrm = set()
for c in df.columns:
    cigars = df[df[c] > scarthreshold].index
    print >> f, '#', c, list(cigars)
    for cxc in it.combinations(cigars[cigars != '720M'], 2):
        s1 = read_csv(path2cigars + '/' + cxc[0], sep = '\t', index_col = 0)
        s2 = read_csv(path2cigars + '/' + cxc[1], sep = '\t', index_col = 0)
        hd = HammingDistance(s1, s2)
        print >> f, cxc, hd
        if hd[0] < hdthreshold:
            if df.loc[cxc[0], c] > df.loc[cxc[1], c]:
                print >> f, '=> '+cxc[1]+' is a sequencing error from '+cxc[0]
                df.loc[cxc[0], c] += df.loc[cxc[1], c]
                df.loc[cxc[1], c] = 0
                scarsrm.add(cxc[1])
            else:
                print >> f, '=> '+cxc[0]+' is a sequencing error from '+cxc[1]
                df.loc[cxc[1], c] += df.loc[cxc[0], c]
                df.loc[cxc[0], c] = 0
                scarsrm.add(cxc[0])

print df.shape, 'shape before cleaning'
df = df[df.sum(axis=1) > minr]
df = 100.*df/df.sum()
print df.shape, 'shape after cleaning'

print scarsrm, 'scars removed in some cells'
for cigar in scarsrm:
    if cigar in df.index:
        print cigar, 'scar removed in a cell, not in others'


df.to_csv(outputfile + '.txt', sep = '\t')
