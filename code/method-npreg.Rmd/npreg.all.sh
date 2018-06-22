#!/bin/bash

# Submit batch job running npreg.R

#SBATCH --time=36:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joycehsiao@uchicago.edu
#SBATCH --partition=broadwl
#SBATCH --account=pi-gilad
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=15
#SBATCH --mem-per-cpu=2GB

set -e

fold=$1

source activate fucci-seq


Rscript npreg.all.R ${fold}
