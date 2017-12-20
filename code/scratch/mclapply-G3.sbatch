#!/bin/bash

# Usage:
#
# sbatch cellcycler-convergence.sbatch ${seed} ${iter}

#SBATCH --time=36:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joycehsiao@uchicago.edu
#SBATCH --partition=broadwl
#SBATCH --account=pi-gilad
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24GB

set -e

seed=$1

source activate fucci-seq

Rscript mclapply-G3.R ${seed}

