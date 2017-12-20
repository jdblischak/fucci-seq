#!/bin/bash

# Usage:
#
# sbatch cellcycler-convergence.sbatch ${seed} ${iter}

#SBATCH --time=36:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joycehsiao@uchicago.edu
#SBATCH --partition=gilad
#SBATCH --account=pi-gilad
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24GB

set -e

seed=$1

iter=$2

source activate fucci-seq

Rscript cellcycleR.R ${seed} ${iter}

