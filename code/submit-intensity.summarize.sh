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

plates=("18511_18855" "18855_19101" "18855_19160" "18870_18511" \
        "18870_18855" "18870_19101" "18870_19160" "19098_18511" \
        "19098_18870" "19098_19160" "19101_18511" "19101_19098" \
        "19160_18870" "18855_19098" "19101_19160" "19160_18511")

#plates=("18511_18855")

dir_intensity="/project2/gilad/fucci-seq/intensities/"
dir_intensity_stats="/project2/gilad/fucci-seq/intensities_stats/"

#dir_proj="/project2/gilad/joycehsiao/fucci-seq/"
#dir_code=${dir_proj}"code/"


for pl in ${plates[@]}; do
#  dir_intensity_pl=${dir_intensity}${pl}/

  sbatch intensity.summarize.sbatch ${dir_intensity} ${dir_intensity_stats} ${pl}

done
