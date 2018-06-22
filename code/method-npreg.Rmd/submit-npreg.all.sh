#!/bin/bash
set -e
#
# Usage:
#
# bash submit-npreg.sh
#

folds=(1 2 3 4 5)

#plates=("18511_18855" "18855_19101")

#dir_intensity="/project2/gilad/fucci-seq/intensities/"
#dir_intensity_stats="/project2/gilad/fucci-seq/intensities_stats/"

for fold in ${folds[@]}; do

  sbatch npreg.all.sh ${fold}

done
