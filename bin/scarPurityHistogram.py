import sys, os
import pandas as pd
import numpy as np
from pandas.io.parsers import read_csv
from collections import Counter
import argparse as argp
import matplotlib.pyplot as plt

try:
    inputdffile = sys.argv[1]
    scartabtxt = sys.argv[2]
    output = sys.argv[3]
except:
    sys.exit('Please, give (1) input norm-filt df; (2) scartab.txt; and (3) output path ')

if not os.path.exists(output):
    os.makedirs(output)

df = read_csv(inputdffile, sep = '\t', index_col = 0)

cigars = list(df.index)

n = {cigar: {i:Counter() for i in range(100)} for cigar in cigars}
with open(scartabtxt) as f:
    for line in f:
        line = line.rstrip().rsplit('\t')
        if line[-1] in cigars:
            cigar = line[-1]
            seq = line[2].rsplit(' ')[0]
            seq = seq[::-1]
            for i in range(len(seq)):
                n[cigar][i].update(seq[i])

n = {cigar: pd.DataFrame(n[cigar]) for cigar in cigars}
n = {cigar: n[cigar].fillna(0) for cigar in cigars}

nt = np.array(['A','C','T','G','N'])
for cigar in n:
    n[cigar] = n[cigar][n[cigar].columns[n[cigar].sum() > 0]]
    n[cigar] = 100.*n[cigar]/n[cigar].sum()
    n[cigar] = n[cigar].transpose()
    n[cigar].index.name = 'loc'
    for nucleotide in nt:
        if nucleotide not in n[cigar].columns:
            n[cigar][nucleotide] = 0
    n[cigar]['seq'] = [n[cigar].columns[n[cigar].loc[i].max()==n[cigar].loc[i]][0] for i in n[cigar].index]
    n[cigar].to_csv(output + '/' + cigar, sep = '\t')

    fig = plt.figure()
    bottom = np.zeros(len(n[cigar]))
    for nt in ['A', 'C', 'G', 'T', 'N']:
        plt.bar(range(len(n[cigar])), n[cigar][nt], bottom = bottom, width = 1)
        bottom += n[cigar][nt]
    art = []
    lgd = plt.legend(['A','C','G','T','N'], loc = 9, bbox_to_anchor = (0.5, -0.15), ncol=5)
    art.append(lgd)
    plt.title(cigar)
    plt.xlabel('position')
    plt.ylabel('nt %')
    fig.savefig(output + '/' + cigar + '.pdf', additional_artist=art, bbox_inches='tight')
sys.exit()

