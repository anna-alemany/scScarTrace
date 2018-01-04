
zcat ${fq}*R1*fastq.gz > ${out}_R1.fastq
zcat ${fq}*R2*fastq.gz > ${out}_R2.fastq

python bin/concatenator.py --fqf ${out} --cbcfile var/bc_scarsc.csv --umifirst --cbchd ${hd} --lenumi 3
cells=384

${path2bwa}/bwa mem -t 8 var/lintrace_histone-GFP_ERCC92.fa ${out}_cbc.fastq > ${out}.sam

python bin/tablator.py --samin ${out} --out ${out} --celnum 384 
