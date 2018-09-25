#!/bin/bash
set -e

inds=("NA19098" "NA18511" "NA18870" "NA19101" \
        "NA18855" "NA19160")

for ind in ${inds[@]}; do
#  dir_images_processed=${dir_scratch}"images_processed"/${pl}/
  sbatch trendfilter.sbatch ${ind}
done

