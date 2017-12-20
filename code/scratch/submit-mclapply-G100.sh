#!/bin/bash
set -e
#
# Usage:
#   bash submit-cellcycler-convergence.sh
#

# iter=({1..20})
# seed=(2  3  5  7 11 13 17 19 23 29 31 37 41 43 47 53 59 61 67 71)

seed=(2)

for ii_seed in ${seed[@]}; do

sbatch mclapply-G100.sbatch ${ii_seed}

done
