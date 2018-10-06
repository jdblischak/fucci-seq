#!/bin/bash

# Submit geo.R to prepare FASTQ files for GEO submission

dir="/project2/gilad/fucci-seq/fastq-combined"
repo="/home/jdblischak/fucci-seq"

for fq in $dir/*/*fastq.gz
do
 name=`basename $fq`
 echo $name
 sbatch -p broadwl --job-name=$name --mem=1g --wrap="Rscript $repo/code/geo.R $fq"
 sleep 1s
done
