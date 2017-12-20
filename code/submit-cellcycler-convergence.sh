#!/bin/bash
set -e
#
# Usage:
#   bash submit-cellcycler-convergence.sh
#

# iter=({1..20})
# seed=(2  3  5  7 11 13 17 19 23 29 31 37 41 43 47 53 59 61 67 71)

iter=(1)
seed=(2)

for i in ${iter[@]}; do

  ii_seed=${seed[i]}
  sbatch cellcycler-convergence.sbatch ${ii_seed} $i

done

