#!/bin/bash

# Usage:
#   sbatch fucciMask.matrixToStats.sbatch
#

#SBATCH --time=36:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joycehsiao@uchicago.edu
#SBATCH --partition=gilad
#SBATCH --account=pi-gilad
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=5GB

set -e

plate=$1

source activate fucci-seq

Rscript fucciMask.matrixToStats.R
