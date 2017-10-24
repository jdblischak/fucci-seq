#!/bin/bash

#SBATCH --time=36:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joycehsiao@uchicago.edu
#SBATCH --mem=15G
#SBATCH --partition=broadwl

set -e

dir_images_data_pl=$1
dir_code=$2
dir_images_processed=$3
dir_intensities=$4
pl=$5

source activate fucci-seq

Rscript ${dir_code}/fucciMask.wrapper.R ${dir_images_data_pl} ${dir_code} ${dir_images_processed} ${dir_intensities} ${pl}
