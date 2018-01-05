# map fastq files with scar reads to reference GFP sequence
# The script requires 3 input parameters: (1) root for fastq R1 and R2 name; (2) root for output file name; (3) path to bwa.
# The script assumes a pair end scar library, where Read 1 contains a 3-nucleotide long UMI (useless in principle), followed by an 8-nucleotide long cell-specific barcode.
# Read 2 contains scar information. By default, no Hamming distance correction is performed on cell-specific barcoes (parameter cbchd)
# The script assumes that each library is made of 384 different cells

# check input parameters
if [ $# -ne 3 ]
then
  echo "Please, give the following inputs in this order"
  echo "1) fastq R1,R2 root"
  echo "2) root for output file name"
  echo "3) path 2 bwa"
fi

fq=$1
out=$2
path2bwa=$3

# merge lanes for each library, when necessary, and unzip fastq.gz files
zcat ${fq}*R1*fastq.gz > ${out}_R1.fastq
zcat ${fq}*R2*fastq.gz > ${out}_R2.fastq

# Select reads with proper cell barcode and produce new ${out}_cbc.fastq file with cell-ID information for each read
python bin/concatenator.py --fqf ${out} --cbcfile var/bc_scarsc.csv --umifirst --cbchd 0 --lenumi 3
cells=384

# Map ${out}_cbc.fastq with bwa
${path2bwa}/bwa mem -t 8 var/lintrace_histone-GFP_ERCC92.fa ${out}_cbc.fastq > ${out}.sam
