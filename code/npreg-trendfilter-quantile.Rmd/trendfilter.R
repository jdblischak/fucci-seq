#!/usr/bin/env Rscript

# Description:
#   This script was written for extracting summary information
#  from the intensity values.
#
# To run this script on the command line
# sbatch

args <- commandArgs(trailingOnly = TRUE)

for (i in seq_along(args)) {
  cat(sprintf("Argument %d: %s", i, args[i]), sep = "\n")
}


library(parallel)
library(genlasso)

dir <- "/project2/gilad/joycehsiao/fucci-seq"

log2cpm.quant <- readRDS(file.path(dir, "output/npreg-trendfilter-quantile.Rmd/log2cpm.quant.rds"))

source(file.path(dir,"code/utility.R"))

#fit.trend <- mclapply(1:nrow(log2cpm.quant), function(g) {
fit.trend <- mclapply(1:100, function(g) {
  fit.trendfilter.generic(log2cpm.quant[g,],
                              polyorder = 2)
}, mc.cores=25)

saveRDS(fit.trend,
        file.path(dir,"output/npreg-trendfilter-quantile.Rmd/fit.trend.quant.rds"))


