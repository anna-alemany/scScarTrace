import sys, os
import argparse as argp
import gzip
from collections import Counter
from Scar import *
import pickle
import pandas as pd

parser = argp.ArgumentParser('GFP Reads per cell')
parser.add_argument('--sam', help = 'input sam.gz file')
parser.add_argument('--out', help = 'root output file')
args = parser.parse_args()

samfile = args.sam
outfile = args.out

def hasRvsPrimer(read):
    primer = 'GGCCCCGTGCTGCTGCCCGAC'
    seq = read[9][-len(primer):]
    hd = 0
    for i in range(len(primer)):
        if primer[i] != seq[i]:
            hd += 1
    return hd <= 3


pile = {}
with gzip.open(samfile) as f:
    for line in f:
        if line[0] == '@':
            continue
        read = line.rstrip().rsplit('\t')
        if read[1] == '16' and 'GFP' in read[2] and hasRvsPrimer(read): 
            cell = int(read[0].rsplit(':')[-1])
            seq = read[9]
            cigar = read[5]
            pos = read[3]
            s = scar(cigar, seq, pos)
            try:
                pile[s].update([cell])
            except:
                pile[s] = Counter([cell])


pickle.dump(pile, open(outfile + '.pickle', 'w'))

df = pd.DataFrame(pile)
df = df.fillna(0)
df = df.transpose()
df.to_csv(outfile + '.tsv', sep = '\t')


