#!/bin/bash

# Submit batch job running run_methods.eva.Rs

#SBATCH --time=36:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joycehsiao@uchicago.edu
#SBATCH --partition=broadwl
#SBATCH --account=pi-mstephens
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=2GB

set -e

ncores=$1
ngenes=$2

source activate fucci-seq

Rscript job_run_methods.train.random.R ${ncores} ${ngenes}
