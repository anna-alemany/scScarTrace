import sys, os
import pandas as pd
import numpy as np
from pandas.io.parsers import read_csv
from collections import Counter
import argparse as argp

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
    n[cigar]['seq'] = [n[cigar].columns[n[cigar].loc[i].max()==n[cigar].loc[i]][0] for i in n[cigar].index]
    n[cigar].to_csv(output + '/' + cigar, sep = '\t')

sys.exit()

