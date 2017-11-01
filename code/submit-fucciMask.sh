#!/bin/bash
set -e
# The following commands capture information about the command line
# arguments:
#
# $1, $2, etc. - the first, second, etc., argument
# $* - all the arguments
# $# - the total number of arguments
#
# Usage:
#
# ./script.sh ex1 ex2 ex3
#
#plates=("18511_18855" "18855_19101" "18855_19160" )
#plates=("18870_18511" \
#        "18870_18855" "18870_19101" "18870_19160" "19098_18511" \
#        "19098_18870" "19098_19160" "19101_18511" "19101_19098" \
#        "19160_18870")

# all plates for which wells do not have duplicated images
#plates=("18511_18855" "18855_19101" "18855_19160" "18870_18511" \
#        "18870_18855" "18870_19101" "18870_19160" "19098_18511" \
#        "19098_18870" "19098_19160" "19101_18511" "19101_19098" \
#        "19160_18870")

#plates=("19101_19160" "19160_18511")
#plates=("18855_19098")

# plates for which some wells have duplicated images need to be analyzed separately
# plates ("18855_19098" "19101_19160" "19160_18511" )

dir_images_data="/project2/gilad/fucci-seq/images_curated/"
dir_intensities="/project2/gilad/fucci-seq/intensities/"

dir_proj="/project2/gilad/joycehsiao/fucci-seq/"
dir_code=${dir_proj}/"code/"

dir_scratch="/scratch/midway2/joycehsiao/fucci-seq/"

for pl in ${plates[@]}; do
  dir_images_data_pl=${dir_images_data}/${pl}/
  dir_images_processed=${dir_scratch}/"images_processed"/${pl}/

  sbatch fucciMask.sbatch ${dir_images_data_pl} ${dir_code} ${dir_images_processed} ${dir_intensities} ${pl}
done
