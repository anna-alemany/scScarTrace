import sys, os
import argparse as argp
import gzip
from Bio import pairwise2
import pandas as pd
import re
from Scar import *
from collections import Counter
import pickle
import itertools as it
import multiprocessing

parser = argp.ArgumentParser('GFP Reads per cell')
parser.add_argument('--picklein', help = 'input sam.gz file')
parser.add_argument('--out', help = 'output file')
parser.add_argument('--th', type = int, help = 'threads', default = 8)
args = parser.parse_args()

pile = pickle.load(open(args.picklein))

def getScarFromAlignment(a):
    ref = a[0]
    seq = a[1]
    j = 0
    cigar = []
    
    if ref.count('-') > 0:
        if ref.index('-') == 0:
            j = min( [ref.index(nt) for nt in ['A', 'C', 'T', 'G'] ])
            ref = ref[j:]
            seq = seq[j:]

    i0 = 0
    n = j
    hamming = 0
    for i in range(len(seq)):
        if i0 == 0:
            if ref[i] != '-' and seq[i] == '-':
                n += 1
            elif ref[i] != '-' and seq[i] != '-':
                n += 1
                i0 = 'm'
        
        elif i0 == 'm':
            if ref[i] != '-' and seq[i] != '-':
                if ref[i] != seq[i]:
                    hamming += 1
                n += 1
            elif ref[i] == '-' and seq[i] != '-':
                cigar.append([n, 'M'])
                i0 = 'i'
                n = 1
            elif ref[i] != '-' and seq[i] == '-':
                cigar.append([n, 'M'])
                i0 = 'd'
                n = 1
    
        elif i0 == 'd':
            if ref[i] != '-' and seq[i] == '-':
                n += 1
            elif ref[i] != '-' and seq[i] != '-':
                cigar.append([n, 'D'])
                i0 = 'm'
                n = 1
            else:
                sys.exit('Err at i0=d')

        elif i0 == 'i':
            if ref[i] == '-' and seq[i] != '-':
                n += 1
            elif ref[i] != '-' and seq[i] != '-':
                cigar.append([n, 'I'])
                n = 1
                i0 = 'm'
            else:
                sys.exit('Err at i0=i')

        else:
            sys.exit('Err at scar determination')
    try:
        cigar[-1][0] += n
    except:
        cigar.append([n, 'M'])

    cigar = ''.join([''.join([str(cigar[i][0]), cigar[i][1]]) for i in range(len(cigar)) ])
    return cigar, hamming


def align(s):
    gfp = 'ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA'
    if s.cigar == '76M' or s.cigar == '75M' or s.cigar == '74M':
        cigar = '720M'
        a = 'WT'
    else:
        a = pairwise2.align.globalms(gfp, s.seq, 1, 0.25, -1, -0.1, penalize_end_gaps = False)[0]
        cigar, mm = getScarFromAlignment(a)
    return (s, cigar, a)

print 'pile read!'

pile2 = {}

pool = multiprocessing.Pool(8) #threads


f = open(args.out + '.txt', 'w')
for idx, (sc, cigar, a) in enumerate(pool.imap_unordered(align, pile.keys())):
    print >> f, sc,'\t',cigar
    print >> f, a
    try:
        pile2[cigar].update(pile[sc])
    except:
        pile2[cigar] = pile[sc]
f.close()

pickle.dump(pile2, open(args.out+'.pickle', 'w'))
df = pd.DataFrame(pile2).transpose().fillna(0).astype(int)
df.index.name = 'CIGAR'
df.to_csv(args.out + '.tsv', sep = '\t')

sys.exit()

#globalms(sequenceA, sequenceB, match, mismatch, open, extend) -> alignments

#match is the score to given to identical characters.  
#mismatch is the score given to non-identical ones.
#open and extend are the gap penalties when a gap is opened and extended.  They should be negative.

#alignments is a list of tuples (seqA, seqB, score, begin, end).
#seqA and seqB are strings showing the alignment between the
#sequences.  score is the score of the alignment.  begin and end
#are indexes into seqA and seqB that indicate the where the
#alignment occurs.




