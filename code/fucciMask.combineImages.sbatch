#!/bin/bash

# Usage:
# dir_code="/project2/gilad/joycehsiao/fucci-seq/code/"
# dir_images_processed="/scratch/midway2/joycehsiao/fucci-seq/images_processed/18511_18855/"
# dir_images_combined="/project2/gilad/fucci-seq/images_combined/"
# dir_intensities="/project2/gilad/fucci-seq/intensities/"
# pl="18511_18855"
#
# sbatch fucciMask.combineImages.sbatch ${dir_code} ${dir_intensities} ${dir_images_processed} ${dir_images_combined} ${pl}

#SBATCH --time=36:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joycehsiao@uchicago.edu
#SBATCH --mem=15G
#SBATCH --partition=broadwl

set -e

dir_code=$1
dir_images_wellId=$2
dir_images_processed=$3
dir_images_combined=$4
dir_intensities=$5
pl=$6

source activate fucci-seq

Rscript ${dir_code}/fucciMask.combineImages.R ${dir_images_wellId} ${dir_intensities} ${dir_images_processed} ${dir_images_combined} ${pl}
