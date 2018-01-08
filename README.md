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
 
## Filtering 
1. *python bin/findThresholds-QCplots.py outfile_rootname.scartab.tsv* <br/>
 This script takes as an input the file *outfile_rootname.scartab.tsv* and provides possible thresholds to filter out cells where scScarTrace did not work and noisy scars that arise due to sequencing and mapping artifacts. The script provides several output files (examples for all these files can be found in the folder "examples"):<br/>
  * ![sumXcell.txt](https://github.com/anna-alemany/scScarTrace/tree/master/examples/sumXcell.txt): text file with the total number of reads detected per cell. If a cell is absent from the list this implies that no read at all was detected. <br/>
  * ![sumXcell.pdf](https://github.com/anna-alemany/scScarTrace/tree/master/examples/sumXcell.pdf): barplot of the total number of reads detected per cell.<br/>
  * ![sumXcell.hst](https://github.com/anna-alemany/scScarTrace/tree/master/examples/sumXcell.hst): normalized histogram of the total number of reads detected per cell.<br/>
  * ![sumXcell.fit](https://github.com/anna-alemany/scScarTrace/tree/master/examples/sumXcell.fit): text file with results from fitting the normalized histogram of the total number of reads per cell with a double Gaussian function. <br/>
  * ![sumXcell-histo.pdf](https://github.com/anna-alemany/scScarTrace/tree/master/examples/sumXcell-histo.pdf): plot showing the histogram of the total number of detected reads per cell and the correspoinding fit to a double Gaussian function. The x-axis is shown in log10 scale. The plot gives a value _x_ for the minimum of the double Gaussian function, which can be used as a threshold to filter out cells with less reads than _10^x_.<br/>
  * ![sumXscar.txt](https://github.com/anna-alemany/scScarTrace/tree/master/examples/sumXscar.txt): text file with number of reads for each detected scar, identified with the cigar code. <br/>
  * ![sumXscar.hst](https://github.com/anna-alemany/scScarTrace/tree/master/examples/sumXscar.hst): text file with the histogram of the number of reads per scar. 
  * ![sumXscar.fit](https://github.com/anna-alemany/scScarTrace/tree/master/examples/sumXscar.fit): text file with results from fitting the histogram of total number of reads per scar to a three domain piecewise-defined linear function. The third domain corresponds to scars that are highly observed and therefore are real scars; the second domain contains scars whose sequence is a few nucleotides away from scars that are present in the third domain, hence these scars are due to sequencing errors of read scars; the first domain corresponds to scars that are hardly observed and therefore are generated by sequencing noise. The border between the second and the third domain _y_ provides a reasonable threshold to filter out scars that are seen with less that _10^y_ reads.<br/>
  * ![sumXscar.pdf](https://github.com/anna-alemany/scScarTrace/tree/master/examples/sumXscar.pdf): plot showing the sorted number of reads for each scar. 
  * ![sumXscar-histo.pdf](https://github.com/anna-alemany/scScarTrace/tree/master/examples/sumXscar-histo.pdf): plot showing the normalized histogram of number of reads per scar and the fit to the three domain piecewise-defined linear function. This profile is general for all scar libraries. 
 
2. *python bin/filter-normalize.py outfile_rootname.scartab.tsv 10^x 10^y outfile_rname* <br/>
 This script takes as an input the scartab.tsv file, followed by the read threshold for cells, the read threshold for scars, and a root name for the output file. As a guide, we recommend to use _~10^x_ as a threshold to filter out cells with less reads, and _~10^y_ as a threshold to filter out scars with less reads. However, the user can explore different thresholds. As an output, the script generates two files:
 * outfile_rname_filter.txt: text file with a table of raw reads for each scar per cell.
 * outfile_rname_norm-filter.txt: text file with a table of normalized reads for each scar per cell. Normalized reads for each cell have been obtained by dividing the number of reads per scar by the total number of reads per cell, and multiplying the result by 100 to have a percentage. 

3. *python bin/scarPurityHistogram.py outfile_rname_norm-filter.txt outfile_rootname.scartab.txt outscardir_name*
 The script takes as an input the outfile_rname_norm-filter.txt file and the outfile_rootname.scartab.txt in order to provide an histogram of sequence content per each cigar. Hence, we can check how many reads for each nucleotide we have in each position of the read. Therefore, when in a given position there are two nucleotides observed with a 50-50 frequency, we know that we have to be careful with this cigar since it codes two different scars. All histograms are stored in the *outscardir_name* directory (which the script creates in case it does not exist). 

4. *python bin/cleanScarErrors.py outfile_rname_norm-filter.txt outscardir_name outfile2_rname*
 

## Clone extraction
2. Cluster
3. Clean noisy scars
4. Final clustering
5. Copy number of each scar
