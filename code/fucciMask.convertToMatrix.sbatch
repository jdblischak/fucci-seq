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
#SBATCH --partition=broadwl
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=20GB
#SBATCH --nodes=3

set -e

plate=$1

source activate fucci-seq

Rscript fucciMask.convertToMatrix.R ${plate}
