#!/bin/bash
set -eux

# sbatch script to test speed of verifyBamID using different number of CPUs.
#
# Usage:
#
# 1 CPU:
# sbatch --tasks-per-node=1 --mem=16G --job-name=test-verify-16G-1cpu --output=test-verify-16G-1cpu.out --partition=gilad /home/jdblischak/fucci-seq/code/test-verifybamid-speeds.sbatch
#
# 4 CPU:
# sbatch --tasks-per-node=4 --mem=16G --job-name=test-verify-16G-4cpu --output=test-verify-16G-4cpu.out --partition=gilad /home/jdblischak/fucci-seq/code/test-verifybamid-speeds.sbatch
#
# 8 CPU:
# sbatch --tasks-per-node=8 --mem=16G --job-name=test-verify-16G-8cpu --output=test-verify-16G-8cpu.out --partition=gilad /home/jdblischak/fucci-seq/code/test-verifybamid-speeds.sbatch

source activate fucci-seq

time=`date +%Y-%m-%d-%H-%M-%S`

echo $time

verifyBamID --vcf /home/jdblischak/fucci-seq/data/snps.grch37.exons.vcf.gz \
            --bam  /project2/gilad/fucci-seq/fucci-bam-verify/20170905/20170905-B02-sort.bam \
            --best --ignoreRG \
            --out /home/jdblischak/fucci-seq/20170905-B02-$time
