#!/bin/bash

# Submit batch job running intensity.summarize.R for one plate

#SBATCH --time=36:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joycehsiao@uchicago.edu
#SBATCH --mem=15G
#SBATCH --partition=broadwl

set -e

dir_intensity=$1
dir_intensity_stats=$2
plate=$3

source activate fucci-seq

Rscript intensity.summarize.R ${dir_intensity} ${dir_intensity_stats} ${plate}
