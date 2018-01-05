# scScarTrace

Abstract for the package. 
In what follows, I will describe each script, its function and its use. Now it is only an item list but will get better

## Initial files
Initially, we have either 8 fastq.gz files for each pair end library, that read as follows:
  * fastqfile_rootname_L001_R1.fastq.gz
  * fastqfile_rootname_L002_R1.fastq.gz
  * fastqfile_rootname_L003_R1.fastq.gz
  * fastqfile_rootname_L004_R1.fastq.gz
  * fastqfile_rootname_L001_R2.fastq.gz
  * fastqfile_rootname_L002_R2.fastq.gz
  * fastqfile_rootname_L003_R2.fastq.gz
  * fastqfile_rootname_L004_R2.fastq.gz
  
In case files from different lanes have already been merged, then we have 2 fastq files:
  * fastqfile_rootname_R1.fastq.gz
  * fastqfile_rootname_R2.fastq.gz

## Mapping
1. *mapFastq.sh fastqfile_rootname outfile_rootname path2bwa* <br/>
  This script takes three input parameters:<br/>
  a) fastqfile_rootname: common part of the name of all fastq.gz files to map <br/>
  b) outfile_rootname: desired name for the outptu files that will be generated <br/>
  c) path2bwa: path to bwa software <br/>
  The script unzips R1 and R2 fastq files and merges files from different lanes. Next, it produces a new fastq file (cbc.fastq) with reads with proper cell-specific barcodes. Finally, it maps the cbc.fastq file using bwa. 
 
2. *gzip outfile_rootname.sam* <br/>
 In order to proceed, it is required to zip the sam file.

3. *python bin/readSAMpileup.py --sam outfile_rootname.sam.gz --out outfile_rootname.pileup* <br/>
 This script goes through all the reads in the sam file, and selects those who have been mapped to the GFP in the proper strand and contain the primer used in the nested PCR for their amplification (GGCCCCGTGCTGCTGCCCGAC) with a Hamming distance equal or less than 3. Then, it counts how many times a given read/scar is seen in each cell. As an output if produces a tabular separated file (tsv) and a pickle file containing the table. 

4. *python bin/realignScars.py --picklein outfile_rootname.pileup.pickle --out outfile_rootname.scartab --th 8* <br/>
 The script remaps pileup reads using the biopython function pairwise2.align.globalms with match, mismatch, open and extend parameters equal to 1, 0.25, -1 and -0.1, respetively. This allows to correct for sequencing errors and pool together reads that belong to the same scar for the same cell. From now on, scars are defined using cigar codes, and the sequence is not used any more. As outptu files we get: <br/>
 a) outfile_rootname.scartab.txt: lists all pileuped reads and corresponding mapping and assigned cigar code <br/>
 b) outfile_rootname.scartab.pickle: pickle version of the previous list <br/>
 c) outfile_rootname.scartab.tsv: table containing number of reads per cell for each scar, denoted now as a cigar <br/>
 

## Clone extraction
1. Filter
2. Cluster
3. Clean noisy scars
4. Final clustering
5. Copy number of each scar
