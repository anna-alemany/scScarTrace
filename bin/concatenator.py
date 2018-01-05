# Reads R1.fastq and R2.fastq files, selects reads with proper cell barcode and produces a new _cbc.fastq file.

import sys, os
import itertools as it
import argparse as argp
import numpy as np

#### function to identify cells from barcodes, allowing some edit distances ####
def getCell(bc, bc2cell, hdmax):
    try:
        cell = bc2sample[bc]
    except:
        cell = 'None'
    return cell
    
def expandBCset(d, hdmax):
    if hdmax == 0:
        nt = ['N']
        hdmax = 1
    else:
        nt = ['N', 'C', 'T', 'G', 'A']

    dg = {}
    for seq in d:
        try: 
            dg[seq].append(d[seq])
        except:
            dg[seq] = [d[seq]]
        i = 0
        while i < hdmax:
            i += 1
            comb = [''.join(l) for l in it.product(nt, repeat = i)]
            for c in comb:
                for p in it.permutations(range(len(seq)), i):
                    s0 = seq
                    for j in range(i):
                        s0 = s0[:p[j]] + c[j] + s0[p[j]+1:]
                    try:
                        dg[s0].append(d[seq])
                    except:
                        dg[s0] = [d[seq]]
    dg2 = dg.copy()
    for seq in dg:
        if len(set(dg[seq])) > 1:
            del dg2[seq]
        else:
            dg2[seq] = dg[seq][0]

    return dg2


#### check input variables ####
parser = argp.ArgumentParser(description = 'Concatenates bcread to bioread qname.')
parser.add_argument('--fqf', help = 'Fastq files names, without _Rx.fastq')
parser.add_argument('--bcread', '-bcr', help = 'read where to find the barcode (umi+cell)', choices = ['R1', 'R2'], default = 'R1')
parser.add_argument('--bioread', '-bior', help = 'read where to find biological information', choices = ['R1', 'R2'], default = 'R2')
parser.add_argument('--lencbc', '-lcbc', help = 'cell barcode length (integer)', type = int, default = 8)
parser.add_argument('--lenumi', '-lumi', help = 'umi length (integer)', type = int, default = 6)
parser.add_argument('--umifirst', help = 'logical variable: umi before cel barcode', action = 'store_true')
parser.add_argument('--cbcfile', '-cbcf', help = 'cell specific barcode file. Please, provide full name')
parser.add_argument('--cbchd', help = 'collapse cell barcodes with the given hamming distance', type = int, default = 0)
args = parser.parse_args()

fqr = args.fqf
bcread = args.bcread
bioread = args.bioread
lcbc = args.lencbc
lumi = args.lenumi
umifirst = args.umifirst
cbcfile = args.cbcfile
hd = args.cbchd

#### Define input fastq files ####
fq1 = fqr + '_R1.fastq'
fq2 = fqr + '_R2.fastq'

if not os.path.isfile(fq1) or not os.path.isfile(fq2):
    print 'fastq files not found'
    sys.exit()

#### Read barcodes ####
d = {}
bcf = open(cbcfile)
for line in bcf.read().rsplit('\n'):
    line = line.rstrip().rsplit('\t')
    if len(line) == 2:
        d[line[0]] = line[1]
        if len(line[0]) != lcbc:
            print 'ERROR: barcode length provided does not match dataset'
            sys.exit()
bcf.close()

bc2sample = expandBCset(d, hd)

#### Do the job ####

fout = open(fqr + '_cbc.fastq', 'w+')
nt = 0
ns = 0
with open(fq1) as f1, open(fq2) as f2: 
    for idx, (l1, l2) in enumerate(it.izip(f1, f2)):
        l1, l2 = l1.rstrip().rsplit(' ')[0], l2.rstrip().rsplit(' ')[0]
        l = np.mod(idx,4)
        if l == 0:
            n1, n2 = l1, l2
            if not n1 == n2:
                sys.exit('fastq files not syncrhonized (@name)')
        if l == 1:
            s1, s2 = l1, l2
        if l == 2:
            p1, p2 = l1[0], l2[0]
            if not p1 == p2 == '+':
                print l1, l2
                print p1, p2
                sys.exit('fastq files not synchronized (+)')
        if l == 3:
            q1, q2 = l1, l2
            if len(q1) != len(s1) or len(q2) != len(s2):
                sys.exit('phred and read length not mathch!')

            if bcread == 'R1':
                bcseq = s1[:lumi+lcbc]
                s1 = s1[lumi+lcbc:]
                q1 = q1[lumi+lcbc:]
            elif bcread == 'R2':
                bcseq = s2[:lumi+lcbc]
                s2 = s2[lumi+lcbc:]
                q2 = q2[lumi+lcbc:]
            if not umifirst:
                celbc = bcseq[:lcbc]
                umi = bcseq[lcbc:]
            else:
                celbc = bcseq[lumi:]
                umi = bcseq[:lumi]
            cell = getCell(celbc, bc2sample, hd)
            if cell == 'None':
                continue
            ns += 1
            if bioread == 'R1':
                print >> fout,  '\n'.join([':'.join([n1, bcseq, umi, celbc, cell]), s1, p1, q1])
            elif bioread == 'R2':
                print >> fout, '\n'.join([':'.join([n2, bcseq, umi, celbc, cell]), s2, p2, q2])
nt = (idx+1)/4
fout.close()

#### LOG ####
fout = open( fqr + '.log', 'w')
print >> fout, '=> to generate cbc file <='
print >> fout, 'fastq file:', fqr
print >> fout, 'full barcode in:', bcread
print >> fout, 'biological read in:', bioread
print >> fout, 'cell specific barcode length:', lcbc
print >> fout, 'umi length:', lumi
print >> fout, 'umi goes first:', umifirst
print >> fout, 'total sequenced reads:', nt
print >> fout, 'reads with proper barcodes:', ns, 1.0*ns/nt
fout.close()
