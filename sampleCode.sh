# ScarTrace: scar extraction and clone identification
# download data from GEO, with GSE102990 accession number. Store it in a data-folder.

libname=data/CB-20170330C-11502-WKM-R2-010-Scar
outname=wkm-r2-010

# Scar extraction

## map scar library with bwa
gunzip data/${libname}.fastq.gz

chmod +x ./bin/mapFastq.sh
./bin/mapFastq.sh data/$libname /usr/local/bin/

gzip data/{libname}.fastq
rm data/{libname}_cbc.fastq

## raw scar extraction
python bin/readSAMpileup.up --sam data/${libname}.sam.gz --out ${outname}.pileup
python bin/realignScars.py --picklein ${outname}.pickle --out ${outname}.scartab --th 8

# Scar filtering

## Quality-check plots
python bin/findThresholds-QCplots.py ${outname}.scartab.tsv

## Read filtering and normalization
python bin/filter-normalize.py ${outname}.scartab.tsv 765 3165 ${outname}-c765s3165 # we filter out cells with less than 765 total reads and scars with less than 3165 reads

## Scar cleaning
python bin/scarPurityHistogram.py ${outname}-c765s3165-norm-filter.txt ${outname}.scartab.txt scardir
python bin/cleanScarErrors.py ${outname}-c765s3165-norm-filter.txt scardir ${outname}-c765s3165-clean-norm-filter

# Clone extraction

## Merge several scar tables
## python bin/mergeDF.py 2 ${outname1}-cxsy-clean-norm-filter lib1  ${outname1}-cxsy-clean-norm-filter lib2 outer merged_file_name

## Hierarchical clustering and further filtering
python bin/hierarchicalClustering.py ${outname}-c765s3165-clean-norm-filter.txt 10 ${outname}-c765s3165-cnf-hcl10 hcl y
python bin/cleanCellScarfraction.py ${outname}-c765s3165-cnf-hcl10_df.txt 3.5 ${outname}-c765s3165-cnf-hcl10-f35 y

## Scar-presence based clustering
python bin/HDclustering.py ${outname}-c765s3165-cnf-hcl10-f35_df.txt ${outname}-c765s3165-cnf35f y

# Clone merging

## automatic scar pooling
python bin/automatPool_HDClustering.py ${outname}-c765s3165-cnf35f_HD.txt ${outname}-c765s3165-cnf35f_HD-apool y

## remove small clones and set wt to clone 0
python bin/rmSamllCl_HDpool.py ${outname}-c765s3165-cnf35f_HD-apool.txt ${outname}-c765s3165-cnf35f_HD-apool_rm y



